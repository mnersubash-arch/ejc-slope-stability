"""Focused reference and boundary tests. This executes Python, never VBA."""
import math, json, csv
from pathlib import Path
import EJC_Eq421_port as e
from Thesis_Columns_port import Config, fastphi, solve

def main():
    tests=[]
    c=Config()
    ph=e.ejc_phi(400,c.sci,c.mb,c.s,c.a,0,1)
    ck=e.ejc_check(400,ph,c.sci,c.mb,c.s,c.a,0,1)
    for key in ("s1","s3","angular","ejc","closure"):
        for value in (math.nan,math.inf,None):
            bad=dict(ck);bad[key]=value
            got=e.ejc_accept(bad)
            tests.append(dict(test=f"{key}={value}",pass_test=got!="ADMISSIBLE",status=got))
        bad=dict(ck);bad.pop(key)
        tests.append(dict(test=f"missing {key}",pass_test=e.ejc_accept(bad)!="ADMISSIBLE"))
    for kw in [dict(tol_ang=math.nan),dict(tol_ejc=math.inf),dict(tol_ang=True),dict(tol_ang=1),dict(tol_ejc=1e300)]:
        tests.append(dict(test=str(kw),pass_test=e.ejc_accept(ck,**kw)=="REJECTED: INPUT"))
    for flag in (math.nan,2,None):
        tests.append(dict(test=f"ordering {flag}",pass_test=e.ejc_accept(dict(ck,ordering=flag))=="REJECTED: ORDERING"))
    # Warm solver versus the scanned solver, including negative normal stresses.
    for a in (.5,.509923):
        c=Config(a=a)
        for n in (-97.9,-92.,-80.,-40.,0.,20.,100.,400.,800.,5000.):
            fast=fastphi(n,c)
            full=None; root_used=None
            for ix in range(1,len(e.brackets(n/1000,c.sci,c.mb,c.s,c.a,0,1))+1):
                candidate=e.ejc_phi(n,c.sci,c.mb,c.s,c.a,0,1,tol=1e-12,root_index=ix)
                if isinstance(candidate,float) and e.ejc_accept(e.ejc_check(n,candidate,c.sci,c.mb,c.s,c.a,0,1))=="ADMISSIBLE":
                    full=candidate;root_used=ix;break
            ok=isinstance(fast,float) and isinstance(full,float) and abs(math.radians(fast-full))<=1e-10
            tests.append(dict(test=f"HB fast/scanned admissible tangent a={a} n={n}",pass_test=ok,fast_deg=fast,scan_deg=full,observed_candidate_index=root_used))
    for a in (.5,.509923):
        r,_=solve(Config(a=a,waterAngle=50,waterHeight=360),3,80)
        ref=1.2053339231344706 if a==.5 else 1.1634915212969383
        tests.append(dict(test=f"primary-port wet 80x80 a={a}",pass_test=r["status"]=="ACCEPTED" and abs(r["F"]-ref)<=1e-8*ref,result=r,reference=ref))
    out=Path(__file__).resolve().parents[1]/"evidence/extended_port_checks.json"
    out.write_text(json.dumps(tests,indent=2))
    print(sum(t["pass_test"] for t in tests),"of",len(tests),"checks pass")
    raise SystemExit(0 if all(t["pass_test"] for t in tests) else 1)

if __name__=="__main__":main()
