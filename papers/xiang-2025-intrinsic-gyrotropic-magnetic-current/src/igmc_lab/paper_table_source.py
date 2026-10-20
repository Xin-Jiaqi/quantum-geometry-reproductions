"""Independent, visually checked transcription of SM pp. 8--10 Tables I/II.

This is source-data transcription, NOT inferred from our numerical MPG bases.
Every nonzero source cell is a linear combination of signed parameters.
Source uses monoclinic unique axis y; source table is stored in source xyz.
The toolkit snapshot instead uses unique axis z for monoclinic 2/m families;
call compare_table with the explicitly declared y<->z permutation.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable
import numpy as np

@dataclass(frozen=True)
class PaperRow:
    id: str
    channel: str
    sm_page: int
    symbols: tuple[str, ...]
    pattern: str

# The arrays below are directly transcribed from the printed SM tables;
# primes use ASCII apostrophe and overbars use leading minus.
# Different written order in Table III does not affect membership.
C_SOURCE_ROWS = (
    PaperRow('C01','C',8,('1',"11'"),'general'),
    PaperRow('C02','C',8,('2',"21'","2'"),'unique_y_2'),
    PaperRow('C03','C',8,('m',"m1'","m'"),'unique_y_m'),
    PaperRow('C04','C',8,('222',"2221'","2'2'2"),'diag3'),
    PaperRow('C05','C',8,('mm2',"mm21'","m'm2'","m'm'2"),'off_xy'),
    PaperRow('C06','C',8,('4',"41'","4'",'3',"31'",'6',"61'","6'"),'uniaxial_gyro'),
    PaperRow('C07','C',8,('-4',"-41'","-4'"),'uniaxial_roto'),
    PaperRow('C08','C',8,('422',"4221'","4'22'","42'2'",'32',"321'","32'",'622',"6221'","6'22'","62'2'"),'uniaxial_diag'),
    PaperRow('C09','C',8,('4mm',"4mm1'","4'm'm","4m'm'",'3m',"3m1'","3m'",'6mm',"6mm1'","6'mm'","6m'm'"),'uniaxial_hall'),
    PaperRow('C10','C',8,('-42m',"-42m1'","-4'2'm","-4'2m'","-42'm'"),'roto_diag'),
    PaperRow('C11','C',9,('23',"231'",'432',"4321'","4'32'"),'isotropic'),
)
D_SOURCE_ROWS = (
    PaperRow('D01','D',9,('1',"-1'"),'general'),
    PaperRow('D02','D',9,('2',"m'","2/m'"),'unique_y_2'),
    PaperRow('D03','D',9,("2'",'m',"2'/m"),'unique_y_m'),
    PaperRow('D04','D',9,('222',"m'm'2","m'm'm'"),'diag3'),
    PaperRow('D05','D',9,("2'2'2",'mm2',"-42'm'"),'off_xy'),
    PaperRow('D06','D',9,("m'm2'","m'mm"),'off_yz'),
    PaperRow('D07','D',9,('4',"-4'","4/m'",'3',"-3'",'6',"-6'","6/m'"),'uniaxial_gyro'),
    PaperRow('D08','D',9,("4'",'-4',"4'/m'"),'uniaxial_roto'),
    PaperRow('D09','D',9,('422',"4m'm'","-4'2m'","4/m'm'm'",'32',"3m'","-3'm'",'622',"6m'm'","-6'm'2","6/m'm'm'"),'uniaxial_diag'),
    PaperRow('D10','D',10,("4'22'","4'm'm",'-42m',"4'/m'm'm"),'roto_diag'),
    PaperRow('D11','D',10,("42'2'",'4mm',"-4'2'm","4/m'mm","32'",'3m',"-3'm","62'2'",'6mm',"-6'm2'","6/m'mm"),'uniaxial_hall'),
    PaperRow('D12','D',10,('23',"m'-3'",'432',"-4'3m'","m'-3'm'"),'isotropic'),
)
ALL_SOURCE_ROWS = C_SOURCE_ROWS+D_SOURCE_ROWS

# Each pattern is represented by explicit source matrix basis. Entries use
# the source's independent symbolic free components, not solver pivots.
def matrix_basis(pattern: str) -> np.ndarray:
    eye = np.eye(3)
    e=lambda i,j: np.outer(eye[i],eye[j])
    patterns={
      'general': [e(i,j) for i in range(3) for j in range(3)],
      'unique_y_2': [e(0,0),e(0,2),e(1,1),e(2,0),e(2,2)],
      'unique_y_m': [e(0,1),e(1,0),e(1,2),e(2,1)],
      'diag3': [e(0,0),e(1,1),e(2,2)],
      'off_xy': [e(0,1),e(1,0)],
      'off_yz': [e(1,2),e(2,1)],
      'uniaxial_gyro': [e(0,0)+e(1,1),e(0,1)-e(1,0),e(2,2)],
      'uniaxial_roto': [e(0,0)-e(1,1),e(0,1)+e(1,0)],
      'uniaxial_diag': [e(0,0)+e(1,1),e(2,2)],
      'uniaxial_hall': [e(0,1)-e(1,0)],
      'roto_diag': [e(0,0)-e(1,1)],
      'isotropic': [e(0,0)+e(1,1)+e(2,2)],
    }
    if pattern not in patterns:raise KeyError(pattern)
    return np.column_stack([m.ravel() for m in patterns[pattern]])

# Toolkit convention unique z vs SM unique y. Physical coordinate change,
# not an antiunitary sign, and only for source mono-axial table templates.
Y_TO_Z = np.array([[1,0,0],[0,0,1],[0,1,0]],dtype=float)
def source_in_toolkit_setting(row: PaperRow) -> np.ndarray:
    basis=matrix_basis(row.pattern)
    if row.pattern.startswith('unique_y_'):
        return np.column_stack([ (Y_TO_Z@b.reshape(3,3)@Y_TO_Z.T).ravel()
                                for b in basis.T ])
    return basis

def row_lookup(rows=ALL_SOURCE_ROWS):
    lookup={}
    for row in rows:
        for symbol in row.symbols:
            key=(row.channel,symbol)
            if key in lookup:
                raise ValueError(f'duplicate SM table MPG: {key} {row.id} {lookup[key].id}')
            lookup[key]=row
    return lookup


def compare_source_table(groups) -> tuple[list[dict],dict]:
    """Full equality of projector onto symbolic tensor space, not just dim."""
    from .toolkit_bridge import tensor_basis
    paper=row_lookup()
    output=[]
    allsymbols={g.symbol for g in groups}
    for channel in ('C','D'):
        extra=sorted(k[1] for k in paper if k[0]==channel and k[1] not in allsymbols)
        if extra:raise AssertionError(f'SM names missing from registry {channel}: {extra}')
        for group in groups:
            key=channel,group.symbol
            source_row=paper.get(key)
            pred=tensor_basis(group,channel)
            if source_row is None:
                target=np.empty((9,0))
                rowid=None
                page=None
            else:
                target=source_in_toolkit_setting(source_row)
                rowid=source_row.id
                page=source_row.sm_page
            # compare two orthogonal projectors, invariant under free-parameter relabeling
            def projector(b):
                if b.shape[1]==0:return np.zeros((9,9))
                q,_=np.linalg.qr(b)
                return q@q.T
            distance=float(np.max(np.abs(projector(target)-projector(pred))))
            actual=int(pred.shape[1]);expected=int(target.shape[1])
            output.append(dict(mpg_number=group.number,symbol=group.symbol,channel=channel,
                        sm_table_row=rowid,sm_page=page,sm_dimension=expected,
                        independently_calculated_dimension=actual,
                        max_projector_residual=distance,
                        status=('PASS' if actual==expected and distance<1e-8 else 'FAIL')))
    result={
      'source_rows':len(ALL_SOURCE_ROWS),
      'source_covered_groups_C':sum(len(r.symbols) for r in C_SOURCE_ROWS),
      'source_covered_groups_D':sum(len(r.symbols) for r in D_SOURCE_ROWS),
      'all_mpg_response_pairs':len(output),
      'matched_mpg_response_pairs':sum(r['status']=='PASS' for r in output),
      'mismatches':[r for r in output if r['status']!='PASS'],
      'max_projector_residual':max(x['max_projector_residual'] for x in output),
    }
    return output,result
