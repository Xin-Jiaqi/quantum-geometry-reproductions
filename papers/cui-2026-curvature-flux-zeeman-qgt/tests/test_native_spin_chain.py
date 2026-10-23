"""Round34 official native Wannier90 formatted spn, U, Udis, V, eig + HR tests."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
from flux_lab.native_spin_chain import (parse_spn_formatted,parse_u,parse_u_dis,
      parse_v,parse_eig,reconstruct_spin,check_against_native_hr)

def fixture(path,*,disentangle=True,bloch_phase=True):
    root=Path(path);root.mkdir(parents=True,exist_ok=True)
    nk=3;nw=2;nb=4 if disentangle else 2
    coords=np.array([[k/3,0.,0.] for k in range(nk)])
    sx=np.array([[0,1],[1,0]],complex)
    sy=np.array([[0,-1j],[1j,0]],complex)
    sz=np.diag([1.,-1.]).astype(complex)
    Sb=np.array([np.kron(np.eye(nb//2),p) for p in (sx,sy,sz)])
    rng=np.random.default_rng(20261023)
    q,_=np.linalg.qr(rng.normal(size=(nb,nw))+1j*rng.normal(size=(nb,nw)))
    E=np.array([[-.41+.04*np.cos(2*np.pi*i/3),
                 -.10+.03*np.sin(2*np.pi*i/3),
                 .26+.03*np.cos(2*np.pi*i/3),
                 .62+.02*np.sin(2*np.pi*i/3)][:nb]
               for i in range(nk)])
    UU=[];UD=[];VV=[];S=[];H=[]
    for ik in range(nk):
        theta=.24+.13*ik
        u=np.cos(theta)*np.eye(2)-1j*np.sin(theta)*sy
        d=np.diag(np.exp(1j*(ik+1)*np.linspace(.13,.73,nb))) if bloch_phase else np.eye(nb)
        Udis=d.conj().T@q
        if not disentangle:
            # Isolated bands have no Udis file: store the FULL Bloch gauge in U.
            u=Udis@u
            V=u
        else:
            V=Udis@u
        Sb_k=np.array([d.conj().T@s@d for s in Sb])
        UU.append(u);UD.append(Udis);VV.append(V);S.append(Sb_k)
        H.append(V.conj().T@np.diag(E[ik])@V)
    UU=np.array(UU);UD=np.array(UD);VV=np.array(VV);S=np.array(S);H=np.array(H)
    paths={n:root/("seed"+ext) for n,ext in
           (("spn",".spn.fmt"),("u","_u.mat"),("dis","_u_dis.mat"),
            ("v","_v.mat"),("eig",".eig"),("hr","_hr.dat"),("r","_r.dat"))}
    spn=["TEST GENERATED, NO REAL DFT","%d %d"%(nb,nk)]
    for ik in range(nk):
        for m in range(nb):
            for n in range(m+1):
                for a in range(3):
                    z=S[ik,a,n,m]
                    spn.append(f"{z.real:.16e} {z.imag:.16e}")
    paths["spn"].write_text("\n".join(spn)+"\n")
    def write_mat(name,M):
        nr,cols=M.shape[1:];lines=["GENERATED, NO REAL MATERIAL",f"{nk} {cols} {nr}"]
        for ik in range(nk):
            lines.append("")
            lines.append(" ".join(f"{k:.16e}" for k in coords[ik]))
            lines.extend(f"{z.real:.16e} {z.imag:.16e}"
                         for z in M[ik].reshape(-1,order="F"))
        paths[name].write_text("\n".join(lines)+"\n")
    write_mat("u",UU)
    if disentangle:
        write_mat("dis",UD)
        write_mat("v",VV)
    paths["eig"].write_text("\n".join(f"{m+1} {ik+1} {E[ik,m]:.16e}"
                         for ik in range(nk) for m in range(nb))+"\n")
    # Independent exact 3-point Fourier: H_R=(1/3) sum_k H(k)e^{-i2pi kR}.
    Rs=[(-1,0,0),(0,0,0),(1,0,0)];deg=[2,1,2]
    Hr=np.array([sum(H[ik]*np.exp(-2j*np.pi*coords[ik,0]*R[0])
                 for ik in range(nk))/nk for R in Rs])
    pos=np.zeros((3,3,nw,nw),complex)
    pos[1,0]=np.diag([.11,.32]);pos[1,1]=np.diag([.04,.20])
    hr=["GENERATED hr","2","3"," ".join(str(x) for x in deg)]
    rr=["GENERATED r","2","3"]
    for ir,R in enumerate(Rs):
        for m in range(nw):
            for n in range(nw):
                z=Hr[ir,m,n]*deg[ir]
                hr.append(f"{R[0]} {R[1]} {R[2]} {m+1} {n+1} {z.real:.16e} {z.imag:.16e}")
                v=pos[ir,:,m,n]*deg[ir]
                rr.append(f"{R[0]} {R[1]} {R[2]} {m+1} {n+1} "+
                          " ".join(f"{x:.16e}" for a in v for x in (a.real,a.imag)))
    paths["hr"].write_text("\n".join(hr)+"\n")
    paths["r"].write_text("\n".join(rr)+"\n")
    # Independent operator result, intentionally not using parser.
    Sw=np.array([[VV[ik].conj().T@S[ik,a]@VV[ik]
                  for a in range(3)] for ik in range(nk)])
    return dict(paths=paths,E=E,U=UU,Udis=UD,V=VV,Sb=S,Sw=Sw,H=H,
                nk=nk,nb=nb,nw=nw,coords=coords)

def load(f):
    p=f["paths"]
    spn=parse_spn_formatted(p["spn"]);u=parse_u(p["u"])
    dis=parse_u_dis(p["dis"]) if p["dis"].exists() else None
    v=parse_v(p["v"]) if p["v"].exists() else None
    eig=parse_eig(p["eig"],nk=f["nk"],nband=f["nb"])
    r=reconstruct_spin(spn,u,u_dis=dis,v=v,eig=eig)
    return spn,u,dis,v,eig,r

def flags():
    return dict(use_ws_distance=False,write_ndegen_applied=False,
                wsvec_present=False,transl_inv_full=False)

class NativeSpnWannierChain(unittest.TestCase):
    def test_actual_spn_triangle_component_index_and_matrix_hermiticity(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td);s=parse_spn_formatted(f["paths"]["spn"])
            np.testing.assert_allclose(s["spin_bloch"],f["Sb"],atol=5e-14)
            self.assertEqual(s["kind"],"BLOCH_PAULI_NOT_WANNIER")
            self.assertEqual(s["nband"],4)
            self.assertEqual(s["nk"],3)

    def test_native_column_major_u_udis_v_matrix_parsing(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td);_,u,dis,v,_,_=load(f)
            np.testing.assert_allclose(u["matrix"],f["U"],atol=5e-14)
            np.testing.assert_allclose(dis["matrix"],f["Udis"],atol=5e-14)
            np.testing.assert_allclose(v["matrix"],f["V"],atol=5e-14)
            np.testing.assert_allclose(dis["matrix"]@u["matrix"],v["matrix"],atol=5e-14)

    def test_independently_projected_physical_pauli_spin(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td);_,_,_,_,_,r=load(f)
            np.testing.assert_allclose(r["spin_wannier_k"],f["Sw"],atol=8e-14)
            self.assertTrue(r["physical_pauli_not_orbital_pseudospin"])
            self.assertFalse(r["physical_source_provenance_authenticated"])
            self.assertFalse(r["material_KME_ready"])
            self.assertGreater(np.max(abs(r["spin_wannier_k"][:,0,0,1])),.03)

    def test_independent_spin_and_native_hr_bloch_eigenvalue_oracle(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td);_,_,_,_,_,r=load(f)
            np.testing.assert_allclose(r["H_wannier_k_from_eig"],f["H"],atol=6e-14)
            gate=check_against_native_hr(r,f["paths"]["hr"],f["paths"]["r"],flags=flags())
            self.assertLess(gate["h_from_eig_vs_hr_max_eV"],1e-12)
            self.assertFalse(gate["physical_spin_source_provenance_authenticated"])
            self.assertEqual(gate["original_paper"],"PARTIAL")

    def test_bloch_gauge_phase_covariance_with_matching_u_and_spn(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(Path(td)/"gauge",bloch_phase=True)
            g=fixture(Path(td)/"plain",bloch_phase=False)
            s=load(f)[-1];t=load(g)[-1]
            np.testing.assert_allclose(s["spin_wannier_k"],t["spin_wannier_k"],atol=8e-14)
            self.assertGreater(np.max(abs(f["Sb"]-g["Sb"])),.08)

    def test_stale_bloch_spn_gauge_fails_phys_spin_even_if_dimensions_match(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(Path(td)/"gauge",bloch_phase=True)
            g=fixture(Path(td)/"plain",bloch_phase=False)
            fresh=load(f)[-1]
            wrong=reconstruct_spin(parse_spn_formatted(g["paths"]["spn"]),
                  parse_u(f["paths"]["u"]),u_dis=parse_u_dis(f["paths"]["dis"]))
            self.assertGreater(np.max(abs(fresh["spin_wannier_k"]-wrong["spin_wannier_k"])),.06)

    def test_omitting_disentanglement_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td);s,u,_,_,_,_=load(f)
            with self.assertRaisesRegex(ValueError,"disentanglement"):
                reconstruct_spin(s,u)

    def test_isolated_bands_do_not_need_u_dis(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td,disentangle=False)
            sp,u,dis,v,eig,r=load(f)
            self.assertIsNone(dis);self.assertIsNone(v)
            np.testing.assert_allclose(r["spin_wannier_k"],f["Sw"],atol=5e-14)
            gate=check_against_native_hr(r,f["paths"]["hr"],f["paths"]["r"],flags=flags())
            self.assertLess(gate["h_from_eig_vs_hr_max_eV"],1e-12)

    def test_h_only_wrong_gauge_hr_eigenvalues_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td);_,_,_,_,_,r=load(f)
            p=f["paths"]["hr"]
            lines=p.read_text().splitlines();columns=lines[8].split()
            columns[5]=str(float(columns[5])+.6)
            lines[8]=" ".join(columns);p.write_text("\n".join(lines)+"\n")
            with self.assertRaisesRegex(ValueError,"HR and Bloch"):
                check_against_native_hr(r,p,f["paths"]["r"],flags=flags(),tol_eV=1e-5)

    def test_bad_v_mismatch_detected_even_with_valid_u_and_ud(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td);p=f["paths"]["v"]
            lines=p.read_text().splitlines();ix=4
            parts=lines[ix].split();parts[0]=str(float(parts[0])+.19)
            lines[ix]=" ".join(parts);p.write_text("\n".join(lines)+"\n")
            with self.assertRaises(ValueError):load(f)

    def test_native_u_k_order_mismatch_with_ud_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td);p=f["paths"]["dis"];lines=p.read_text().splitlines()
            lines[3]="0.2700000000 0 0";p.write_text("\n".join(lines)+"\n")
            with self.assertRaisesRegex(ValueError,"k order"):load(f)

    def test_truncated_spn_and_wrong_triangular_count_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td);p=f["paths"]["spn"]
            lines=p.read_text().splitlines()
            p.write_text("\n".join(lines[:-1])+"\n")
            with self.assertRaisesRegex(ValueError,"upper-triangle"):
                parse_spn_formatted(p)

    def test_nonphysical_spin_operator_norm_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td);p=f["paths"]["spn"]
            lines=p.read_text().splitlines()
            lines[2]="2.0000 0.0000";p.write_text("\n".join(lines)+"\n")
            with self.assertRaisesRegex(ValueError,"exceeds one"):parse_spn_formatted(p)

    def test_nonhermitian_diagonal_imaginary_spin_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td);p=f["paths"]["spn"]
            lines=p.read_text().splitlines()
            lines[2]="0.2 0.007";p.write_text("\n".join(lines)+"\n")
            with self.assertRaisesRegex(ValueError,"diagonal"):parse_spn_formatted(p)

    def test_invalid_eig_order_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td);p=f["paths"]["eig"]
            lines=p.read_text().splitlines();lines[0]="2 1 -0.2"
            p.write_text("\n".join(lines)+"\n")
            with self.assertRaisesRegex(ValueError,"ordering"):
                parse_eig(p,nk=3,nband=4)

    def test_native_u_matrix_col_major_negative_control(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td);p=f["paths"]["u"]
            lines=p.read_text().splitlines()
            # Source values are F-order, C-order produces U.T (different for general U).
            m=load(f)[1]["matrix"][0]
            misread=np.asarray([complex(*map(float,s.split())) for s in lines[4:8]]).reshape((2,2))
            self.assertGreater(np.max(abs(m-misread)),.1)

    def test_full_source_authentication_remains_blocked(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td);sp,u,dis,v,eig,r=load(f)
            self.assertEqual(sp["kind"],"BLOCH_PAULI_NOT_WANNIER")
            self.assertEqual(len(sp["source_sha256"]),64)
            self.assertEqual(len(dis["source_sha256"]),64)
            self.assertTrue(r["mathematical_gauge_transform_verified"])
            self.assertFalse(r["physical_source_provenance_authenticated"])
            self.assertFalse(r["material_KME_ready"])
            self.assertEqual(r["original_cui_paper"],"PARTIAL")

if __name__=="__main__":unittest.main()
