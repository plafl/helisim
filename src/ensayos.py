from matematicas import *
from Lynx import modelo
modelo.init()

from pylab import *
rc('text', usetex=True)

tr1 = modelo.trimado(30, 0, 0, 0, 1000, 1.225)
tr2 = modelo.trimado(60, 0, 0, 0, 1000, 1.225)

modelo.hG = 0
modelo.viento_vel = 0
modelo.viento_dir = 0
modelo.P = 99461.0575
modelo.T = 283
modelo.motor.nmotores = 2

modelo.t = 0.0
modelo.x_i = modelo.y_i = 0.0
modelo.z_i = 1000

modelo.trimado(30, 0, 0, 0, 1000, 1.225)
#s = modelo.series(5., 
#        escalon(1., None, tr1['et0p'],  tr2['et0p']),
#        escalon(1., None, tr1['et1sp'], tr2['et1sp']),
#        escalon(1., None, tr1['et1cp'], tr2['et1cp']),
#        escalon(1., None, tr1['etpp'],  tr2['etpp']),
#        0.03, 
#        'u', 'v', 'w', 'p', 'q', 'r', 'th', 'fi', 'ch', 'th0', 'th1s', 'th1c', 'th0T')

s = modelo.series(15., 
        tr1['et0p'],
        escalon(1., None, tr1['et1sp'], tr1['et1sp'] - 0.1),
        tr1['et1cp'],
        tr1['etpp'],
        0.03, 
        'u', 'v', 'w', 'p', 'q', 'r', 'th', 'fi', 'ch', 'th0', 'th1s', 'th1c', 'th0T')

plot(s['t'], s['et1sp'])
xlabel(r'$t$')
ylabel(r'$\eta_{1sp}$')
show()

plot(s['t'], s['th1s'])
xlabel(r'$t$')
ylabel(r'$\theta_{1s}$')
show()

plot(s['t'], s['th'])
xlabel(r'$t$')
ylabel(r'$\theta$')
show()

plot(s['t'], s['fi'])
xlabel(r'$t$')
ylabel(r'$\phi$')
show()

plot(s['t'], s['u'])
xlabel(r'$t$')
ylabel(r'$u$')
show()

plot(s['t'], s['v'])
xlabel(r'$t$')
ylabel(r'$v$')
show()

plot(s['t'], s['w'])
xlabel(r'$t$')
ylabel(r'$w$')
show()

plot(s['t'], s['p'])
xlabel(r'$t$')
ylabel(r'$p$')
show()

plot(s['t'], s['q'])
xlabel(r'$t$')
ylabel(r'$q$')
show()

plot(s['t'], s['r'])
xlabel(r'$t$')
ylabel(r'$r$')
show()

modelo.t = 0.0
modelo.x_i = modelo.y_i = 0.0
modelo.z_i = 1000

tr0 = modelo.trimado(0, 0, 0, 0, 1000, 1.225)
s = modelo.series(15., 
        escalon(1, None, tr0['et0p'], tr0['et0p'] + 0.1),
        tr0['et1sp'],
        tr0['et1cp'],
        tr0['etpp'],
        0.03, 
        'u', 'v', 'w', 'p', 'q', 'r', 'th', 'fi', 'ch', 'th0', 'th1s', 'th1c', 'th0T')

plot(s['t'], s['et0p'])
xlabel(r'$t$')
ylabel(r'$\eta_{0p}$')
show()

plot(s['t'], s['th'])
xlabel(r'$t$')
ylabel(r'$\theta$')
show()

plot(s['t'], s['fi'])
xlabel(r'$t$')
ylabel(r'$\phi$')
show()

plot(s['t'], s['u'])
xlabel(r'$t$')
ylabel(r'$u$')
show()

plot(s['t'], s['v'])
xlabel(r'$t$')
ylabel(r'$v$')
show()

plot(s['t'], s['w'])
xlabel(r'$t$')
ylabel(r'$w$')
show()

plot(s['t'], s['p'])
xlabel(r'$t$')
ylabel(r'$p$')
show()

plot(s['t'], s['q'])
xlabel(r'$t$')
ylabel(r'$q$')
show()

plot(s['t'], s['r'])
xlabel(r'$t$')
ylabel(r'$r$')
show()
