using OverconfidentTraders

ot = OverconfidentTraders

par  = ot.Pars()
beta = ot.findbeta(par)
bpar = ot.CstrBPar(beta, par)

a    = 1
ϵ    = 0

z    = a + ϵ/beta
p    = ot.price(z, par, bpar)