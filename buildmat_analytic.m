function diagonal=buildmat_analytic(bb,hx,hy)
a0=bb(1);
a1=bb(2);
a2=bb(3);
a3=bb(4);
a4=bb(5);
a5=bb(6);
diagonal = 0.3215020576e-5 * (0.2697754697e5 * a4 * a0 * hx ^ 2 + 0.644227871e4 * a5 ^ 2 * hy ^ 4 + 0.1348877349e5 * a2 ^ 2 * hy ^ 2 + 0.2339853519e4 * a3 ^ 2 * hx ^ 2 * hy ^ 2 + 0.2697754697e5 * a0 * a5 * hy ^ 2 + 0.77760e5 * a0 ^ 2 + 0.1348877349e5 * a1 ^ 2 * hx ^ 2 + 0.644227871e4 * a4 ^ 2 * hx ^ 4 + 0.467970705e4 * a4 * a5 * hx ^ 2 * hy ^ 2) / hy / hx;

