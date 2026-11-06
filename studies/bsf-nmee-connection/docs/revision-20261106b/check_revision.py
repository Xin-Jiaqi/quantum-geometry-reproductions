from pathlib import Path
import json,csv,re,math
r=Path(__file__).resolve().parent
def run(root):
 d=json.loads((root/'figure_data.json').read_text());ab={(x['t_meV'],x['U_top_minus_bottom_meV']):x for x in d['km']}
 errors=[]
 for b in d['km_ba']:
  a=ab[(b['t_meV'],-b['U_top_minus_bottom_meV'])]
  errors.extend([abs(a['bottom']+b['top']),abs(a['top']+b['bottom']),abs(a['total']+b['total'])])
 assert max(errors)<1e-12
 assert min(v['common_gap_meV'] for v in d['km']+d['km_ba'])>0
 assert ab[(0.,0.)]['bottom']<0 and ab[(20.,0.)]['bottom']>0
 assert ab[(10.,0.)]['bottom']<0 and ab[(10.,2.)]['bottom']>0
 assert all(v['total']>0 for v in d['km'] if v['t_meV']==10)
 z=list(csv.DictReader((root/'decoupled-bias-check.csv').open()));delta=max(float(x['bottom']) for x in z)-min(float(x['bottom']) for x in z);assert delta<1e-12
 for p in root.glob('*.md'):
  s=p.read_text();assert not any(ord(c)<32 and c not in '\n\t' for c in s),p
  assert s.count('$$')%2==0,p
  for link in re.findall(r'\]\(([^)]+)\)',s):
   if '://' not in link: assert (p.parent/link).exists(),(p.name,link)
 s=(root/'正文讨论稿-修订版.md').read_text()
 assert all(x.count('|')==6 for x in s.splitlines() if x.startswith('|'))
 for n in range(1,4): assert '图 '+str(n)+'｜' in s
 sm=(root/'补充讨论与图.md').read_text()
 for n in range(1,6):assert '图 S'+str(n)+'｜' in sm
 result={'mirror_response_max_abs':max(errors),'decoupled_bias_variation':delta,'min_common_gap_meV':min(v['common_gap_meV'] for v in d['km']+d['km_ba']),'main_captions':3,'supplement_captions':5,'table_columns':5,'local_links':'pass','scope':'Data, symmetry relation, insulating condition and manuscript integrity; not new DFT convergence.'}
 (root/'revision_checks.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':run(r)
