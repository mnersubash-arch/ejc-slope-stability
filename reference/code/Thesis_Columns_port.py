"""Python control-flow mirror of the new primary VBA assembly, for pre-native verification.
Not execution of VBA, and not a replacement for Microsoft Excel testing.
Strength uses the identical equation 4.21 and analytical Newton-on-phi derivative.
"""
from dataclasses import dataclass, asdict
import math, json
import EJC_Eq421_port as ejc

@dataclass
class Config:
    xc: float = 355.7423115464064
    zc: float = 644.4192948645371
    radius: float = 725.2669975187362
    height: float = 360.
    slope: float = math.tan(math.radians(50))
    gamma: float = 27.
    uExtra: float = 0.
    waterAngle: float = 0.
    waterHeight: float = 0.
    strength: int = 1
    closure: int = 3
    nu: float = .271
    sci: float = 77.7
    mb: float = 10*math.exp(-58/28)
    s: float = math.exp(-58/9)
    a: float = .5
    cMC: float = 100.
    phiMC: float = 0.
    tolA: float = 1e-10
    tolE: float = 1e-8
    tolF: float = 1e-8
    checks: int = 1
    maxIter: int = 2000

def coeff(c):
    return {1:(-1,1/c.nu),2:(c.nu,c.nu),3:(0,1)}[c.closure]

def geometry_issue(c):
    """Limits of the adopted toe-cutoff, positive-centre, driving-side grid.

    This is a supported-domain guard, not a generalisation of the grid. Allow only
    1E-5 R of toe-intersection roundoff (the published rounded circle is 2.62 mm).
    """
    values=(c.xc,c.zc,c.radius,c.height,c.slope,c.gamma)
    if not all(isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v) for v in values):
        return "UNSUPPORTED GEOMETRY: finite numeric geometry inputs required."
    if c.radius<=c.zc or c.zc<=0 or c.height<=0 or c.slope<=0 or c.gamma<=0:
        return "UNSUPPORTED GEOMETRY: require R > zc > 0 and positive height, slope and unit weight."
    if c.xc<=0 or c.xc>=c.radius:
        return "UNSUPPORTED GEOMETRY: this driving-side toe-cutoff model requires 0 < xc < R."
    x_toe=c.xc-math.sqrt(c.radius**2-c.zc**2)
    if x_toe < -1e-5*c.radius:
        return "UNSUPPORTED GEOMETRY: the adopted grid would omit slope-side mass between its right bound and the toe."
    return ""

def fastphi(sn, c, guess=-1, tol=1e-12):
    p,q = coeff(c)
    if (p,q)!=(0,1):
        return ejc.ejc_phi(sn,c.sci,c.mb,c.s,c.a,p,q,tol)
    args=(sn/1000,c.sci,c.mb,c.s,c.a,p,q)
    pb=ejc.phi_b(*args[:5])
    if pb is None:return "#NUM!"
    lo,hi=pb+1e-9,math.radians(ejc.HI_DEG)
    if lo>=hi:return "#NUM!"
    def gs(x):return ejc.gstate(x,*args)
    lv,hv=gs(lo),gs(hi)
    if not lv or not hv or lv[0]*hv[0]>0:
        return ejc.ejc_phi(sn,c.sci,c.mb,c.s,c.a,p,q,tol)
    glo=lv[0];x=(lo+hi)/2
    if lo<math.radians(guess)<hi:x=math.radians(guess)
    def ar(x):
        check=ejc.ejc_check(sn,math.degrees(x),c.sci,c.mb,c.s,c.a,p,q,tol)
        return check.get("angular") if check.get("domain")==1 else None
    for _ in range(200):
        gv=gs(x)
        if gv is None:break
        g,dg,_=gv;r=ar(x)
        if r is not None and dg!=0 and abs(r)<=tol and abs(g/dg)<=tol:return math.degrees(x)
        if glo*g<=0:hi=x
        else:lo=x;glo=g
        xn=x-g/dg if dg else (lo+hi)/2
        if not lo<xn<hi:xn=(lo+hi)/2
        step=abs(xn-x);x=xn;r=ar(x)
        if r is not None and abs(r)<=tol and step<=tol:return math.degrees(x)
    return ejc.ejc_phi(sn,c.sci,c.mb,c.s,c.a,p,q,tol)

def elements(c, dim, n):
    issue=geometry_issue(c)
    if issue: raise ValueError(issue)
    assert dim in (2,3) and type(n) is int and 2<=n<=20000
    assert dim!=3 or n<=400
    x0=c.xc-math.sqrt(c.radius**2-c.zc**2);x1=c.xc-c.radius
    dx=(x0-x1)/n;dy=c.radius/n if dim==3 else 1.
    out=[]
    for i in range(n):
        x=x0-(i+.5)*dx
        for j in range(n if dim==3 else 1):
            y=(j+.5)*dy if dim==3 else 0.
            q=c.radius**2-(x-c.xc)**2-y*y
            if q<=0 or x>0:continue
            zb=c.zc-math.sqrt(q);h=min(c.height,-x*c.slope)-zb
            if h<=0:continue
            ah=dx*dy;g=math.sqrt(q)/c.radius;rho=math.sqrt(c.radius**2-y*y)
            u=c.uExtra
            if c.waterHeight>0:u+=9.81*max(0,min(c.waterHeight,-x*math.tan(math.radians(c.waterAngle)))-zb)
            out.append(dict(x=x,y=y,zb=zb,h=h,ah=ah,g=g,rho=rho,area=ah/g,kap=(c.xc-x)/rho,
                            weight=c.gamma*h*ah,pressure=u,sn=(c.gamma*h-u)*g*g,phi=-1,
                            label="NOT EVALUATED"))
    return out

def strength(c,e):
    p,q=coeff(c)
    if c.strength==0:e["phi"]=c.phiMC;e["cohesion"]=c.cMC;return True
    ph=fastphi(e["sn"],c,e["phi"],min(1e-12,c.tolA/10) if (p,q)==(0,1) else c.tolA)
    if isinstance(ph,str):e["label"]="NO CONVERGED EQ.4.21 ROOT";return False
    b=c.mb*e["sn"]/(1000*c.sci)+c.s
    if b<=0:e["label"]="OUTSIDE HB DOMAIN";return False
    r=math.radians(ph);ht=(b/(1+math.sin(r)/c.a))**c.a
    e["phi"]=ph;e["cohesion"]=1000*c.sci*ht*math.cos(r)/2-e["sn"]*math.tan(r)
    chk=ejc.ejc_check_state(e["sn"],e["cohesion"],ph,c.sci,c.mb,c.s,c.a,p,q,c.tolA)
    e["label"]=ejc.ejc_accept(chk,c.tolA,c.tolE)
    return e["label"]=="ADMISSIBLE"

def validate(c,els,F,rebuild=True):
    p,q=coeff(c);angles=[];errors=[];equil=[];normals=[];flags=[]
    resistance=driving=momentR=momentD=0.
    for e in els:
        t=math.tan(math.radians(e["phi"]));ma=e["g"]+e["kap"]*t/F
        if ma<=0:return dict(pass_state=False,evaluated=False)
        if rebuild:e["finalN"]=(e["weight"]-e["pressure"]*e["ah"]-e["cohesion"]*e["area"]*e["kap"]/F)/(e["area"]*ma)
        nn=e["finalN"];tau=e["cohesion"]+nn*t
        eq=abs(e["weight"]-((nn+e["pressure"])*e["ah"]+tau*e["area"]*e["kap"]/F))/e["weight"]
        label="ADMISSIBLE"
        if tau<0 or eq>1e-9:label="REJECTED: EQUILIBRIUM/SHEAR"
        if c.strength==1:
            chk=ejc.ejc_check_state(nn,e["cohesion"],e["phi"],c.sci,c.mb,c.s,c.a,p,q,c.tolA)
            lbl=ejc.ejc_accept(chk,c.tolA,c.tolE)
            if lbl!="ADMISSIBLE":label=lbl
            angles.append(abs(chk.get("angular",math.inf)));errors.append(abs(chk.get("ejc",math.inf)))
        e["tau"]=tau;e["label"]=label;e["equilibrium"]=eq
        equil.append(eq);normals.append(nn);flags.append(label=="ADMISSIBLE")
        resistance+=tau*e["area"];driving+=e["weight"]*e["kap"]
        momentR+=tau*e["area"]*e["rho"];momentD+=e["weight"]*e["kap"]*e["rho"]
    globalRes=abs(resistance/driving-F)/F
    return dict(pass_state=all(flags) and globalRes<=c.tolF,evaluated=True,admissible=sum(flags),
                maxAngle=max(angles,default=0),maxEJC=max(errors,default=0),maxEq=max(equil),
                globalRes=globalRes,moment=momentR/momentD,minNormal=min(normals),tension=sum(x<0 for x in normals))

def solve(c,dim,n):
    issue=geometry_issue(c)
    if issue:
        return dict(status="NOT ACCEPTED",F=None,active=0,iterations=0,evaluated=False,reason=issue),[]
    els=elements(c,dim,n)
    if not els:
        return dict(status="NOT ACCEPTED",F=None,active=0,iterations=0,evaluated=False,reason="No active sliding elements."),[]
    limit=-1000*c.s*c.sci/c.mb;F=1.;proj=0
    if c.strength==1:
        for e in els:
            if e["sn"]<=limit:e["sn"]=limit+max(1e-9,.001*abs(limit));proj+=1
    for it in range(1,c.maxIter+1):
        fails=0;drive=0.;resist=0.;move=0.
        for e in els:
            fails+=not strength(c,e);drive+=e["weight"]*e["kap"]
        if fails:return dict(status="NOT ACCEPTED",F=None,active=len(els),iterations=it,reason=f"{fails} failed strengths"),els
        for e in els:
            t=math.tan(math.radians(e["phi"]));ma=e["g"]+e["kap"]*t/F
            if ma<=0:return dict(status="NOT ACCEPTED",F=None,reason="NONPOSITIVE m_alpha"),els
            sf=(e["cohesion"]*e["ah"]+(e["weight"]-e["pressure"]*e["ah"])*t)/ma
            resist+=sf
            raw=(e["weight"]-e["pressure"]*e["ah"]-e["cohesion"]*e["area"]*e["kap"]/F)/(e["area"]*ma)
            if c.strength==1 and raw<=limit:raw=limit+max(1e-9,.001*abs(limit));proj+=1
            move=max(move,abs(raw-e["sn"]));e["sn"]=raw
        fn=resist/drive
        if fn<=0:return dict(status="NOT ACCEPTED",F=None,reason="NONPOSITIVE F"),els
        df=abs(fn-F)/fn
        if df<=min(1e-12,c.tolF/10) and move<=1e-9:
            v=validate(c,els,fn)
            if v["pass_state"]:
                return dict(status="ACCEPTED" if c.checks==1 else "DIAGNOSTIC ONLY",
                            F=fn if c.checks==1 else None,internalF=fn,active=len(els),iterations=it,
                            dF=df,projections=proj,**v),els
        F=fn
    return dict(status="NOT ACCEPTED",F=None,reason="ITERATION LIMIT"),els

BENCHMARKS=[
    ("HB dry .5",{},2.5536684200376,2.31025502757599),
    ("HB wet .5",dict(waterAngle=50.,waterHeight=360.),1.20012184467095,1.03587083245705),
    ("HB dry .509923",dict(a=.509923),2.47603808311863,2.24206326532995),
    ("HB wet .509923",dict(a=.509923,waterAngle=50.,waterHeight=360.),1.15842745018062,1.00046615780691),
    ("Published",dict(xc=207.28,zc=586.53,radius=622.08),2.14193739499583,2.00387341089056),
    ("MC phi0",dict(strength=0),.255354533635107,.190263799071471),
]
if __name__=="__main__":
    results=[]
    for label,kw,f3,f2 in BENCHMARKS:
        for dim,n,expected in [(3,20,f3),(2,40,f2)]:
            cfg=Config(**kw);r,es=solve(cfg,dim,n)
            ok=r["status"]=="ACCEPTED" and abs(r["F"]-expected)<=1e-9*expected
            results.append(dict(test=label,dim=dim,mesh=n,expected=expected,passed=ok,result=r,inputs=asdict(cfg)))
            print(label,dim,r.get("F"),"PASS" if ok else "FAIL",flush=True)
    for closure in (1,2):
        for dim,n in [(3,20),(2,40)]:
            r,_=solve(Config(closure=closure),dim,n)
            results.append(dict(test=f"closure {closure} withheld",dim=dim,passed=r["status"]=="NOT ACCEPTED" and r["F"] is None,result=r))
            print("closure",closure,dim,r,flush=True)
    cfg=Config(a=.509923,waterAngle=50,waterHeight=360)
    r,els=solve(cfg,3,20);els[0]["finalN"]+=1
    mutation=validate(cfg,els,r["F"],rebuild=False)
    results.append(dict(test="actual final stress +1 kPa",passed=not mutation["pass_state"],result=mutation))
    r,_=solve(Config(checks=0),3,20)
    results.append(dict(test="checks off withheld",passed=r["status"]=="DIAGNOSTIC ONLY" and r["F"] is None,result=r))
    unsupported=[
        ("original MC example: unsupported negative centre",dict(xc=-1.45,zc=3.24,radius=4.46,height=2.1,slope=1/3.5,gamma=19,strength=0,cMC=7.2,phiMC=30)),
        ("positive centre: omitted toe-side mass",dict(xc=207.28,zc=300,radius=622.08)),
    ]
    for label,kw in unsupported:
        for dim,n in [(3,20),(2,40)]:
            cfg=Config(**kw);r,els=solve(cfg,dim,n)
            ok=r["status"]=="NOT ACCEPTED" and r["F"] is None and r["active"]==0 and not r["evaluated"] and not els and r["reason"].startswith("UNSUPPORTED GEOMETRY:")
            results.append(dict(test=label,dim=dim,passed=ok,result=r,inputs=asdict(cfg)))
            print(label,dim,r["reason"],"PASS" if ok else "FAIL",flush=True)
    from pathlib import Path
    Path("primary_engine_port_results.json").write_text(json.dumps(results,indent=2))
    raise SystemExit(0 if all(r["passed"] for r in results) else 1)
