from pylab import *

def puttext(x0,y0,name='some name',ce=0,cq=0,size=18.,col='w'):
    text(x0,y0,name,size=14,color=col)
    text(x0,y0-0.1,r'$c_E$ = %d'%ce,size=size,color=col)
    text(x0,y0-0.2,r'$c_Q$ = %d'%cq,size=size,color=col)

cdict = {
'red'  :  ((0., 0.1, 0.1), (10.1/11, 0.5, 0.5), (1., 0.8, 0.8)),
'green':  ((0., 0.1, 0.1), (10.1/11, 0.5, 0.5), (1., 0.8, 0.8)),
'blue' :  ((0., 0.1, 0.1), (10.1/11, 0.5, 0.5), (1., 0.8, 0.8))
}

my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

e = linspace(0,1,100)
q = linspace(0,1,100)

ctrl_e = 0.5*tanh(200.0*(e-0.5)/5.5)+0.5
ctrl_q = 0.5*tanh(200.0*(q-0.5)/5.5)+0.5

s = zeros((100,100))

for i in range(100):
  for j in range(100):
    s[j,i] = 0.1*max(ctrl_e[i],ctrl_q[j]) + 0.9*(1-ctrl_e[i])*ctrl_q[j] - (1-ctrl_q[j])*(1-ctrl_e[i])*10.0

f = figure(figsize=(6,5))

pcolor(e,q,s,cmap=my_cmap)
xlabel('E')
ylabel('Q')
title(r'$v$ / $v_0$')
puttext(0.1,0.35,name='akinetes',ce=0,cq=0)
puttext(0.6,0.35,name='heterocysts',ce=1,cq=0)
puttext(0.1,0.85,name='recruitives',ce=0,cq=1)
puttext(0.6,0.85,name='vegetatives',ce=1,cq=1)
plot([0,1],[0.5,0.5],'k--')
plot([0.5,0.5],[0,1],'k--')
colorbar(ticks=[-10,-5.0,-1.0,0.1,1.0])
savefig('sinking_factor.pdf',dpi=300.)
close()
