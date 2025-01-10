
//    [ H + D_x    0        A_e^T  A_d^T  A_i^T    0     0   |                              ] [ x   ] = [ -f_x ]
//    [    0     D_x          0      0    -I      -I     I   |                              ] [ s   ] = [ -f_s ]
//    [ A_e       0        -D_e      0     0       0     0   |I   -I                        ] [ λ_e ] = [ -g_e ]
//    [ A_d       0          0       0     0       0     0   |                              ] [ λ_d ] = [ -g_d ]
//    [ A_i      -I          0       0     0     -D_i    0   |        I  -I                 ] [ λ_i ] = [ -g_i ]
//    [   0     Zl_i         0       0     0     Sl_i    0   |                              ] [ zl  ] = [ -cl  ]
//    [   0    -Zu_i         0       0     0      0    Su_i  |                              ] [ zu  ] = [ -cu  ]
//    -----------------------------------------------------------------------------------------------------------
//    [                      I                                 D_x           -I             ] [ pe  ] = [ -g_pe ]
//    [                     -I                                     D_x          -I          ] [ ne  ] = [ -g_ne ]
//    [                                    I                           D_x         -I       ] [ pi  ] = [ -g_pi ]
//    [                                   -I                              D_x         -I    ] [ ni  ] = [ -g_ni ]
//    [                                                       Zpe            Spe            ] [ zpe ] = [ -cpe ]
//    [                                                           Zne           Sne         ] [ zne ] = [ -cne ]
//    [                                                               Zpi          Spi      ] [ zpi ] = [ -cpi ]
//    [                                                                  Zni          Sni   ] [ zni ] = [ -cni ]
// 
//  zpe = Spe^{-1} (-cpe - Zpe pe)
//  => lam_e  + (D_x+Spe^{-1}Zpe^{-1}) pe = -g_pe - Spe^{-1}cpe
//  pe = (D_x+Spe^{-1}Zpe^{-1})^{-1} (-g_pe - Spe^{-1}cpe - lam_e)

//    [ H + D_x    0        A_e^T  A_d^T  A_i^T    0     0  ] [ x   ] = [ -f_x ]
//    [    0     D_x          0      0    -I      -I     I  ] [ s   ] = [ -f_s ]
//    [ A_e       0        -D_ee     0     0       0     0  ] [ λ_e ] = [ -gg_e]
//    [ A_d       0          0       0     0       0     0  ] [ λ_d ] = [ -g_d ]
//    [ A_i      -I          0       0     0     -D_ii   0  ] [ λ_i ] = [ -gg_i]
//    [   0     Zl_i         0       0     0     Sl_i    0  ] [ zl  ] = [ -cl  ]
//    [   0    -Zu_i         0       0     0      0    Su_i ] [ zu  ] = [ -cu  ]
//  D_ee = D_x + (D_x+Spe^{-1}Zpe^{-1})^{-1} + (D_x+Spi^{-1}Zpi^{-1})^{-1}
//  gg_e = -g_e + (D_x+Spe^{-1}Zpe^{-1})^{-1} (-g_pe - Spe^{-1}cpe) - (D_x+Spi^{-1}Zpi^{-1})^{-1} (-g_pi - Spi^{-1}cpi)