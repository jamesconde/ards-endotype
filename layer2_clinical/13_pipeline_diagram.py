import sys; sys.path.insert(0,'/home/james/ards-endotypes')
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from shared import paths

INK,ACC,WARN,MUT="#2b2b2b","#1f6feb","#c0392b","#6b7280"
fig,ax=plt.subplots(figsize=(11.5,6.4)); ax.set_xlim(0,100); ax.set_ylim(0,64); ax.axis("off")

def box(x,y,w,h,t,c=INK,fc="white",fs=8.2,bold=False):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.45",
        fc=fc,ec=c,lw=1.4))
    ax.text(x+w/2,y+h/2,t,ha="center",va="center",fontsize=fs,color=INK,
            fontweight="bold" if bold else "normal",linespacing=1.35)

def arrow(x1,y1,x2,y2,c=MUT,style="-|>",lw=1.3,ls="-"):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle=style,color=c,lw=lw,
        mutation_scale=13,linestyle=ls,shrinkA=2,shrinkB=2))

ax.text(50,61.5,"Pipeline: reproduce, then stress-test under real screening conditions",
        ha="center",fontsize=12,fontweight="bold")

# data row
box(1,50,20,7,"Challenge 2019\n40,336 ICU patients\n2 hospital systems",ACC,"#eef4ff",8,True)
box(24,50,17,7,"eICU demo\n2,520 stays\n186 hospitals",ACC,"#eef4ff")
box(44,50,16,7,"GEO x3\n802 / 144 / 62\nsamples",ACC,"#eef4ff")
box(63,50,16,7,"MGH Olink\n383 subjects\n1,472 assays",ACC,"#eef4ff")
ax.text(91,53.5,"6 sources\n3 modalities",ha="center",fontsize=8,color=MUT,style="italic")

# layer 2 spine
box(1,39,20,6.5,"Cohort: hypoxemic\nresp. failure  n=5,687")
box(1,29,20,6.5,"Reference assignment\n(cohort-standardised)\nNOT deployable",MUT)
box(1,19,20,6.5,"Portable classifier\nkappa 0.852\npasses both checks",ACC,"#eef4ff",8.2,True)

box(25,19,19,6.5,"Transfer A->B\n82.4% vs 1.6%\nscoreable",WARN,"#fff0ee",8.2,True)
box(25,29,19,6.5,"Severity check\n+0.017 AUC over\nmodified SOFA")
box(25,39,19,6.5,"Ablation (simulated)\nkappa 0.869 -> 0.461\nat 6 variables")

box(48,29,20,6.5,"eICU replication\nmortality 16.8% vs 36.6%\np<0.0001",ACC,"#eef4ff",8.2,True)
box(48,39,20,6.5,"APACHE IVa benchmark\n0.700 vs our 0.634\n(external model)")
box(48,19,20,6.5,"Refinement loop\nRR 1.62 -> 2.02\nholdout scored once")

box(72,24,26,12,"ENRICHMENT SIMULATION\n\nScreen 12,494  vs  647,654\nto enrol the same trial\n\n52x, from measurement practice",
    WARN,"#fff0ee",9,True)

# layer 1
box(25,7,19,7,"MARS anchor 87.5%\nclusters recover published\nlabels  ARI 0.385")
box(48,7,20,7,"Frozen signature ->\n2 held-out cohorts\nzero refitting")
box(72,7,26,7,"Olink trajectory\nday-0 class -> 28-day acuity\n40.3% switch class by day 7",ACC,"#eef4ff")

for a in [(11,50,11,45.5),(11,39,11,35.5),(11,29,11,25.5),(21,22,25,22),
          (34.5,25.5,34.5,29),(34.5,35.5,34.5,39),(44,22,48,22),
          (58,25.5,58,29),(58,35.5,58,39),(68,30,72,30),(44,32,48,32)]:
    arrow(*a)
arrow(52,50,34.5,14,MUT,"-|>",1.1,"--")   # GEO -> MARS anchor
arrow(71,50,85,14,MUT,"-|>",1.1,"--")     # Olink -> Olink trajectory
arrow(44,10.5,48,10.5)                    # anchor -> frozen transfer
arrow(68,10.5,72,10.5)                    # transfer -> Olink

ax.text(1,15.5,"LAYER 2 — clinical",fontsize=8.5,color=MUT,fontweight="bold")
ax.text(1,3.5,"LAYER 1 — molecular",fontsize=8.5,color=MUT,fontweight="bold")
ax.plot([0,100],[16.5,16.5],color="#e5e7eb",lw=1)
ax.text(99,1,"every model frozen before transfer; holdouts scored once",
        ha="right",fontsize=7.5,color=MUT,style="italic")

fig.savefig(paths.FIGURES/"pipeline_overview.png",dpi=130,bbox_inches="tight",facecolor="white")
print("wrote pipeline_overview.png")
