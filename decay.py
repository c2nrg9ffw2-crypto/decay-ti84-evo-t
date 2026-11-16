# Decay chain solver for the TI-84 Evo
# Type the number of protons (Z) and neutrons (N).
# The program follows the main decay path on the nuclide chart
# until the nucleus is stable, and shows every step.
# Data: IAEA nuclide chart (ground states), made by decay_make.py
try:
    from ti_system import disp_clr, wait_key
except:
    disp_clr = wait_key = None
try:
    import gc
except:
    gc = None

# --- DATA START (made by decay_make.py, do not edit) ---
NMIN = [1, 0, 1, 1, 2, 2, 2, 3, 3, 4, 5, 6, 6, 8, 8, 10, 10, 12, 12,
    12, 14, 17, 16, 19, 18, 19, 19, 21, 20, 24, 24, 28, 27, 30, 29,
    32, 31, 35, 35, 37, 37, 40, 39, 42, 41, 44, 44, 45, 46, 47, 49,
    52, 52, 55, 54, 57, 58, 60, 63, 62, 65, 67, 67, 67, 70, 70, 73,
    73, 75, 75, 79, 78, 79, 82, 83, 84, 85, 87, 87, 91, 90, 95, 96,
    101, 102, 106, 107, 110, 113, 116, 117, 120, 122, 126, 134, 128,
    137, 136, 138, 141, 141, 143, 146, 149, 149, 150, 152, 153, 155,
    157, 157, 161, 165, 165, 170, 172, 174, 176, 176]
CNT = [6, 7, 8, 9, 11, 15, 15, 15, 16, 19, 20, 21, 23, 23, 23, 22,
    23, 24, 24, 26, 25, 24, 27, 26, 29, 30, 32, 31, 35, 31, 32, 29,
    32, 30, 33, 35, 36, 35, 36, 36, 38, 37, 39, 38, 41, 40, 42, 41,
    41, 42, 42, 40, 42, 40, 43, 41, 41, 41, 38, 41, 39, 38, 40, 40,
    39, 40, 38, 39, 38, 38, 37, 40, 40, 40, 41, 41, 43, 42, 44, 41,
    47, 42, 43, 41, 42, 39, 39, 37, 34, 32, 32, 29, 42, 26, 20, 26,
    19, 21, 21, 18, 19, 17, 15, 15, 15, 16, 14, 19, 16, 17, 16, 15,
    10, 13, 7, 4, 5, 2, 1]
SYM = ("n H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr"
    " Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og")
MODES = (
    "buuuuussbnmuussnbnbnnppsshbnbnaeasbbbnbnnpcpssbbbbnbnbnmqccc"
    "ssbbbbbbbubppccssbbbbbbbbnqqgccsssbbbbbbnmuuppccsbbbbbbbbnbn"
    "bqqcccsssbbbbbbbbbbnbuppeccsbbbbbbbbbbbbnbqqccccsssbbbbbbbbb"
    "bnbnbpcccccsbbbbbbbbbbdbbbbbccccccsssbbbbbbbbbbbbbbpcccccsbb"
    "bbbbbbbbbbbbbpcccccsssbsbbbbbbbbbbubppccccsbsbbbbbbbbbbbbbbb"
    "qcccccsesbsbbbbbbbbbbbbtpuppccccsbsbbbbbbbbbbbbbbbpcccccsess"
    "sbsblbbbbbbbbtdppcccccsbbbbbbbbbbbbbbubqcgcccecsssssbbbbbbbb"
    "bbbbdbpccccccelsbbbbbbbbbbbbbbbbcccccccclesssbbbbbbbbbbbbubb"
    "beucccccccecsbbbbbbbbbbbbbbbubbqccccccccsesssbbbbbbbbbbbbbbb"
    "tbbpucccccccecsbbbbbbbbbbbbbbbbbbbqcccccccccscsssbsbbbbbbbbb"
    "bbbbbbbbuppccccccccscsbbbbbbbbbbbbbbbbbbqcccccccccscsssblbbb"
    "bbbbbbbbbbbtpcccccccecsbsbbbbbbbbbbbbbbbbqccccccccecsesssbsb"
    "bbbbbbbbbbbbbpcccccccccecsbbbbbbbbbbbbbbbbbcccecccccecsesssb"
    "sblbbbbbbbbbbbbbpupcccccccccsbsbbbbbbbbbbbbbbbbbuubecccccccc"
    "ccscsesssbsbbbbbbbbbbbbbbbbppcccccccccecsbbbbbbbbbbbbbbbbbbb"
    "bbcccccccccecsesssbbbbbbbbbbbbbbbbbbbbpccccccccccccsbbbbbbbb"
    "bbbbbbbbbbbbbbcccccccccccecsssbsblbbbbbbbbbbbbbbbbbbcccccccc"
    "ccccsbbbbbbbbbbbbbbbbbbbbbbbbgccccccccccsesssssblbbbbbbbbbbb"
    "bbbbttbbpcccceccccccebbbbbbbbbbbbbbbbbbbbbbbbbpcuccccccccscs"
    "ssssbsbbbbbbbbbbbbbbbbbbbbbpcccccccccccecsbbbbbbbbbbbbbbbbbb"
    "bbbbbbbccccccccccecsesssbsbsbbbbbbbbbbbbbbbbbbbbbppccceccccc"
    "ccccsbsbbbbbbbbbbbbbbbbbbbbbbbgcceccccccccscsesssbsblbbbbbbb"
    "bbbbbbbbbbbpccgcccccccccccecsbbbbbbbbbbbbbbbbbbbbbbbbccccccc"
    "ccccecscsssssssbsbsbbbbbbbbbbbbbbbbpcccccccccccccccecsbsbbbb"
    "bbbbbbbbbbdbbbbaaaaccccccccccecscslsssblblbbbbbbbbbbbubbbapc"
    "ccccccccccccccecsbbbbbbbbbbbbbbubudbbaaaacccccccccceclcsesss"
    "ssblblbbbbbbbbbbbdbbppcccccccccccccccccecsbbbbbbbbbbbbbbbbbb"
    "bccccccccccccccecsclesssssbbbbbbbbbbbbbbbbpuuccccccccecccccc"
    "ccecsbbbbbbbbbbbbbbbbbbcccccccccccccecsclesblbbbbbbbbbbbbbbb"
    "bpuucccccccccccccccccsbbbbbbbbbbbbbbbbbbbbcccccccccccccccecs"
    "sassbsblbbbbbbbbbbbbbecccccccccccccccceebbbbbbbbbbbbbbbbbbbc"
    "ccccccccccccccseaaassbsbsbbbbbbbbbbbbbbpppucccccccccccccecec"
    "lcsbbbbbbbbbbbbbbbbucucccccccccccacaeaesssssbsbbbbbbbbbbbbpu"
    "ucccccccucccccccccececsbbbbbbbbbbbbbbbcccccccccccccccacscses"
    "ssssbbbbbbbbbbbbppcccccccccccccccccccececsbbbbbbbbbbdbbucucc"
    "ccccaaccccececscsesssbsbbbbbbbbbbpppccccccaaccccccccccccecsb"
    "bbbbbbbbbbbcucccaaccccccccececsesssssbsbbbbbbbbbpppcauaaaccc"
    "ccccccccccccecsbbbbbbbbbbbbbuuuccaaccccccccccccececaesssssbb"
    "ubbbbbbbppaacaaccccccccccccccccceesbbbbbbbbbbbbbeaaaaccccccc"
    "cccccccececleslsbsbbbbubbbbbbuppacaccccccccccccccccccecsbbbb"
    "bbbubbbbbbaaaaaaaccccccecccccccecleassssbsbbbbbbbbbbbcpaaaac"
    "ccceccccccccccccccecsbsbbbbbbbbbbbbaaaaaaaaaaacccccccccccccc"
    "aesesssbsbbbbbbbbbbppaaaauccccccccccccccccccecsbbbbbbbbbbbbb"
    "aaaaaaaaaaccccccccccccececsesssssbsbbbbbubbbbbbpaaacaccccacc"
    "ccucccccccccecsbsbbbbbbbbbbbbaaaaaaaacccccccccccccceceelesss"
    "bbbbbbbbbbbbapaaaaaaccccccccccccccccclbabbbbbbbbbdddbaaaaaaa"
    "aaaacacccccccccaaaaaaaaaaabbbubbbbbaaaaaaaaaaaccccccccceaaaa"
    "aaaabbbbbbbbbdaaaaacaacaaacacacaeaaaaaaaaababbbbbbbubaaaaaaa"
    "aaaaaaaacaaaaaaaaabbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaaaaababbbb"
    "bbbbaaaaaaaaaaaaaaaaaaaeabbbbbbbbbubaaaaaaaaaaaaaaaaaaaaaaaa"
    "babbbbbbaaaaaaaaaaaaaaaaacecabbbbbbbbaaaaaauaaaaaaaaeaeaaaaa"
    "babbbbbbbbbbbbbubbbaauaaaaaaeaccceceeabbbbbbbaaacececaeaaaba"
    "babbbauuuuuaeueceeecceeababbbbucaaeaeeaeaaaaaaabfbcaucucefac"
    "ceecaubbbuuuafeacacaeaeaaaaabfbfaaaaeecueceeaaabbbffafaaaaca"
    "caeaafaffaaaaeaececceeeaffuafaaaacaafafufcaaaaaaaaffufufffff"
    "afafafaufufaaaaaaaafuufffuffaaafaffffuauaaaauauaauuaaauauuuf"
    "aaaaaaaaauauauffauauauuuaaafauuuaauaaauauuuaufffaauauuuaaafa"
    "uuuaauuuafafaaauuuaaaaauuuafaaaaacaaaaaaaauaaa")
HL = (
    "Sp00002J00000000a3001{2U3v00002.PH3aOX003L002c0000PKOd00NE00"
    "3uXT8N00fVQjNb3XMs2]2w2wPD6F0000NaNXNSNP3qM!00Me372>3fOYQoSz"
    "0000cwPoPBOeONN!NW00M/2J2uNRSp0000QCP+O|OtOYOFNdNUHc2a3,NCRm"
    "Ru000000QwQjPyPmOSO300DB2k2p2-4MRlT%00QgP+P+PmOxOBN8Mz4vMa00"
    "002?3#OWPhQm000000Q+R-O_OfNmNaNVM}MjMkI9Hk003T4MO)QrZC00U)R5"
    "PbOpNkNyN%NXNTN6M)MTI9Hk3iD6OLOXP(Qg000000SkU/PdOtOmOHOLNaNR"
    "N2N600I900HLOMO)PlQDe|00RwSHP!P0NnNwN/NsOLNRNGN2IH00I8I8NjNw"
    "OZOhPmP+000000T/bCQ2PsPEO)OLN[N%NnNaNSHkIRHGNyOlOmP*Rz00W<X1"
    "QhQ?P{PnP0OnOaOVN%NqNYIB00H(NWOXOePcPq000000X#00S1T;QhQTPkPa"
    "OlOVN}N)0000H6HGOeOpPpPf00ef00S?TCRoQ,Q9PvO>O!OiOV00I8000000"
    "D)NVOTOcPKPi00XB00bN00T%aOS5SrQqQPPcO#OcOW0000IrD)00HBHBOdOt"
    "PcSS00ih00U#U:S#SxRsQnQ9PdO+OvOWNkNQ0000HLNgOVOdO(PM00eK0000"
    "00X{00WGtPScQiQfP/O*OWNcNRIrIr00ILOdO_P4U0U100XyW4V;TDRrQhQN"
    "PqO;ORNgNTNS000000I3NiN+ODOgO/aoT>0000000000S9RrQ%PlPdOfOTN;"
    "NiNcNVIrIR00HfOAOWO=O#S;W>YErC00R:RpQ_Q6OhOrOeOSOXN%NoNXNVIR"
    "IRIrNTNbNxN>OlO.U:S|s4X4000000R.SAQqQBPAO-OiOgOYNxNhNe00IRIr"
    "IrI100NqOVObOxOnT2WRfrYC00T/RoPuP/OnP7ONOnOLONN_N#NiNWNa00Ir"
    "00MUNTNcN!N_OaOpUiSa00ZD000000XKfgSBRlQ2PkPHO(O_OeOWN_NqI6NT"
    "N8Ir000000NtO0OVOjOeU,XrY8Xl00ZdT!RpQwOpPcOfO#OfOiOWOBN=NvNl"
    "NkNcNT00MWN1NYNeNvN.OWOfWWV.00d{000000b600T.V[QqQyQgQ1PqPgPK"
    "O/OtOjObOWNxNe0000I4HzN/OBOfPwRoS(T[Sm00U%00S2V{Q!R%Q<QoQ7P+"
    "PgPcP0O+OtOjOWO4Nn00MRNaNmNuOIOdRyRpUrS[00Y5000000TCsTRyV<Qt"
    "RqQfP|PlPfPBO>OrOcOXN-ItHTO1OcOXQ#R!SuUuW4Tv00S!00U(UASWRvQ%"
    "QiP>PsPjPcO_OqOQONNxIsNTI2NyOYOaRlQ!T,SyY8V/00W/000000Ty00U!"
    "TzQoQzQHP)PiPVO.OiOZIL0000HTO0OYORQ:RzSuT9V|V*XuW[00V*V/TzSf"
    "QkQ%QoQiP+PkPUO,Of00000000NT00NnNwOYQ*QwS{R|WtUX00X?000000ei"
    "00Sxu7S#R,Q%QjP_PfO%OfOm0000I6IU00HLHAOAQqRnR-S)T!U*V]SF00Sx"
    "00V.T-S:R(R1R2QlP-PjO<OqOVO1I6I60000000000N0NcNiNuOQQmQwSrR{"
    "U)Tw00V.00eY000000a200TxT;R*Q#QRPiPdOgOWOBN[NxNTM_Is00I1H*N_"
    "QoQ*R;SxS%Q%U7RnX!X900W]k1SxSuR!R4P.P}PrOxOfOcOXN-N+NmNrNd00"
    "0000NfNhOJQKQVR#RwT%S#X1V-00Xf000000XRaKUvT/SPRnQtPbO%P1OmOf"
    "OXO0N,N,NtNaIU00IBN:N,QkQzRmQOSMS_T:U)W4X>00V|XZT_UzSySpP[P%"
    "O=PfP9O)OvOiOfOWOANoNkNf000000I8N/P/P_Q#Q:S*SUU*T#XxW4000000"
    "fW00XetGU+Q!PlQCPnPtPdPcP2OdOaO8N/NsNeNk0000IBN)P(QdQpRpR:St"
    "T)U)b%gu00dQXBU;TwPtQkPeQCP,PfP<PtPaOpOfOWODN-NnNmNXNd0000Ia"
    "00M.PnPwQoQjSVRvUHSu00cf0000000000V|szSuSqRlRkQ*QSPzPbO{OpOe"
    "OXOBN;N!NmNcNZ00I1N.PmQ5QiQ]R*R?T;UAU.WEfwfweWQlStP[R0SxSRQ*"
    "QqP?PRPQOoOmOaOVO9N:NyNkNcNbNc00IaIZ00PcPdQhQLR/RkT8T!00V}00"
    "00000000XF00U5Y_R;R_Q(QhPlPiPGO<OrOfOaOUO1NzNjNfNZNV00I3NjPf"
    "P:QhRmS1SoS/ScW>U.ZJY200Q:V.QzS#QmRnPyQgP!PsPiPZP4O(OoOcOYOa"
    "ILNwNkNgNZNa000000PaPaQVP=RuR*SxS!W7Uk00W?000000f_00U(00S%U/"
    "RqRyQuQhP,PjPSO-OnOcOWNsN:N%NsNpNl000000IENgPiP-QvQ[RvRuSqSs"
    "TvTwXHS(00Ry00QuWjT>UFP/SzR;RmP%PlPcPEO;OpOeObN+OWN;N!NwNpNi"
    "I*NmPaPsQYQlQ]RoS6SOTETC00UQ00Y}000000pu00V?tTT.T6R#Q{QjP?Pl"
    "PcP4O/OwOnOaObN}OFN:N_00NqNkPvP{QkQsRlRsS2SDS:TEU3UAV}Su00Rm"
    "oFQjS}P=RyPvQsPfQ3PvPoPfPbPKO{OoOnOfObOZOVOGN_I*PcPhP%QBQqQ%"
    "RtR(SpSxU3S=00X?00000000000000V*00X[00W*eXT*TFRwR:R2Q.PfPaO:"
    "OtOeOZOY00HZO(PcO_P)QFQmQtRmR0SIR.S:SvT;R/V/Sv00V|00XaZEW:W9"
    "UpU5S_S%R%RxPEPhPSO)OuOOOcOVN,H4IrK.MgPlP-QoQoRuRrSuS9T.TuWV"
    "U*00W]00q:000000UszDTwu]S)W3SrS|QoQnPpPePgOuOe00IW0000NqL9P2"
    "PpPxQ7PlRnPtRwStSyTyT*R/U%WDXZW;00S)gbU#WpT,U/T8URRoQtQ4PnPM"
    "O#OhOYILIWOP00KtNTOOPAPrPrQfQnR5RkR;S9S_S_U.T)n/U+00XC000000"
    "0000WNw,UqvSR;StQ.QjPhPcO/PbOeOaOJIV0000LtKJO?PeP6QPQjQ;RkR!"
    "QqSAQ!T3RqUNR/V,S-W+Wa00Z7fdW;aLS:SiRlQuPhPiPZO]OrOiOaOWOFO0"
    "N)O%O)PdPiP_P]QtQzRuRySqR,T#SsV]T+00W:vaa20000000000TyW;SxSq"
    "QjQhP,PmPPO{OuOlOcOZOXN,Ne0000PsP[QRQlQqRlR0S2S3SrScTFU9U0SF"
    "U.Sod)kp00V:U0TzStQ/QuQ2P*PdPbO]O*OpOjObOVOF00Pb00P%Q1QcQ{Q("
    "R<R.S%SpT_T!W3U,00UpqPX]00X9qVV-Y9R)SsR2R3P>P)PiPePMP8OqOiOd"
    "OUNQ0000PcPxPvP+PsQzQ.RpRqSGSwS(SsTxRoU5R-00U-W<SwULS(SsRxRx"
    "Q3QoP!P,PnPfO(OpOdOYOc00P100PiP=Q8QiQuRqRmSaSrT7S[UBS.W5T.00"
    "00o}0000W/00T#s{SrSqQ#QvQUP[PbPIO,O(OhOq00PaPoPqQ4Q3QjP=Q]S1"
    "RyQfR?QXS!Q.Y7YHa9ZfZCWPV?T:V+R>S3R#Q:QwQgP;PfP8PaO}O%0000O="
    "00PcP)PtQaQfQ?Q=R*R!SuSpTwSd00YFh4kqpk0000a_00V<00S#UtSVS4Qg"
    "QbP;PoPcPePYPG0000MBNY0000O.PfPxQgQhQnPfQ/PnR!QfWUWHX0XUX,aS"
    "s8a500Z/ZYW=U)T2SxQ;QvQgQIP+PnPhPdOf0000PbIBPmP:P}QlQjRmQ-R_"
    "S%XOV/a!W%fYX?n*Y50000000000U,00R:SZRlQ=QfP;P,PuPBO!0000MF00"
    "00IBPgPnPzO_QhPa00QLT!TuU2T_U,U,V[U:WOWOa!bE00XmWeSRSzR)RvQu"
    "QoQNP>PWPcPC0000O_00PQPnP{QWQ1Q%RlR,R?SNSxT-UPfkUy00Ug00X]00"
    "00000000T,W5SDScQ-R1P*I7PeI70000M.MqOz00P6PoPxP}PmQqRmQ)R#Ru"
    "SrT4TCSsSqS;S*T.SuclS-00V*T>R)R|R%Q}QuQAPwPjI7000000IB00PhPw"
    "P/P)QnQtQfQ+R:S4SzSyT,S>V+T?00Tw00Uz000000W(00UaV>RoR+RmI700"
    "I70000JFJSN}O[P6PQPmP+QLPfQMQrRoR/R=SgSjS.S#T%RuV,UcW#X,00X["
    "Z5V{UhS5SuRsRqILI7ILI7P6IBPgPuP+O!PiQvQ-RpRrR}R?SySqTxSoV]Sx"
    "00X80000000000WD00T(TwSVRyI7I7IDI7I7IaNzOCP6PQ00N}O-Q9QgQhQ*"
    "RnRoR=R*SqR!T8S6V-V=WrWcZ0ZJ00j;WbS,U7S8R.RuR4QoI70000000000"
    "00PkPKNdOXPsP{QjQnQ-Q.RsRnSJRuS*R+U*U#Z5V)o_Xk0000000000XIgK"
    "TuU2R.R!00000000MeOWNQN.PJPg00P!QgQjQ!Q(RnRuS0SJS%S?T>TvUzUg"
    "V[T-Z4Uh00X>WMUmT5SpR|Qp00P[ILPm0000OnMNM}OMO!PcPrQ4P>QoQpQ{"
    "RmRyRySHSRS;S=T.T+W|S?s9X?00u+00Xp00V)XkSrS.00000000I7I7I700"
    "L*LoOWOyP6PgPmQ0P-QMQXQkQkRuRySAS4StSsSzRzU.V|XkXB00W8j[U+V)"
    "R)SnQl00P=Q1PuI7I7I7L-MWM)NbO2OgPHPlPyQFQOQoQrQ<RoR/R)S1SGS#"
    "T%U:U%nMX-o_00000000W=00V+ZkSGS=R%00P=Q100I70000J6NQNpOhOuPN"
    "PwP-QVQKQVQSQzQhRnRpS0SuTET>U(U+UzV:W;W:00Xo00U-T,Q|S9QLQ1IL"
    "QqQg00I700LWLaL=MVM_NUN!OTOxPPPpQ4QfQqQqR2Q|R#SGSwTwT*T-W.U!"
    "l6V}00af000000U.00S/U%RzV;QrQf00I70000LZKNNcNgOX0000PfPqQCQP"
    "QjQlQ;QqR?SqSYSdS-S}T?UAU,V/Y0WX00V|W2T4S*QxRkQ.Q#Q.IL00ILIL"
    "K_KyLTM1MWNRNaOXOmPaPqP!QgQZQ!Q]RoRtR+SRSzT5U9T}bfU!00V|0000"
    "000000XN00S2SYR(S{Q*00IL00000000ILM#NYOkOiPbPwPvQAQfQpPyQ{Rm"
    "RxR!00SlS#S;TwT%T;UEUZV*W1W:00ZO00R?R}R*RwRnRpQ!QtQgQcIL00LI"
    "MkMqNzN.O<O-Q4P;QkQvQ-RmRnR.00SqSuS?SWT-TzU:UsdwV?r0gc000000"
    "T?aDS>U!SpS+RyILQpQk00ILNTKxNVNrN=P5Q4QhQ)RlRqR)S2SiSpS+S>T#"
    "T#U#U!W=WXaNemtGWLRvTuT2SzSRRwRqQ%QrILILILIL00KTMPLXMkMZNcNm"
    "OzOyP/P}RoRsS6SrSvT1S?T_T#WxUJZFb8X]O:IKJXLMMTOZPfR*SpILRtSh"
    "IL00ILILILMSOJNiOoOoOyOyP+QBQ;RpR*SPSgS+S/T%T!UFUgUXOqI3IlLG"
    "LaNnPdR2R:RxR0Q_RnILILILILILMML}M.MtN_N_O]PaQBQcQ<RmR%S8ShS("
    "S-T-U)S(NZIHJJKkLyNoP)R2S)W9S(T%R|SPQpRlQh0000L)NVNSN(N[OpO="
    "PiP(QlQkR5Q_R+R*SzQ(M)H/IzKNMMNeQwS0StS#R,R=Q]RzQ,Q_QoQnP_PQ"
    "N6NWNlN:OgOjPePdP;P#QiQiR#PoMSI9JCKRNGNYQxQ,W/W7W=cHS|ZhR=T!"
    "RrR?QzQzNaNcNhOQOIOuOgPOPAQNOcLoHtJ6KDNgN+P=RvT;W-V+aDUNTuRu"
    "SQRuRyQ<RkRmNQMSMaNWNrNmOaOIPcNgLVI3J6K8MSMXO_PaSdS/W]Z5c]d_"
    "V*joS#X0SNS[R}SjILM!M%NXNUOaMnLHHeI,JtMeM!PLPhRsS[U:V.W?ddV,"
    "X3USS(SgScRxT%LwL=MuNWL.Kz00IvJhKJLkO0OmRlSgTEW_WDayeReaiBgi"
    "Wdi=S%U(00Sw00000000000000000000000000LLKQ00ITJIKdMlNpO/RkR="
    "R{T5SuS>WFY{eQfbV>V]TuStRwRsRxPbRpRrSbS<S!UnS)ZFXMa[dUc*a5em"
    "UBhHUzW/V[Mz0000000000PQQ!00RnR+RxSpR/TwT!U#V>bdU*c<UzT)S]S%"
    "00QsQ{00SJ00T+T:X3X9X{aKa9d0cmgbekTvc}SwQqRx00Qr00Ry00R}R{SL"
    "U7UBWLV;cFZ=YET?TC00Sp00PlNbQ-R4R;R.SqSyT1V.T>YEbVa4c1ZCW[Xa"
    "TySrQ1QLQnQqQ+RlSQ00S(T#UlV-Y}W_Y9XGS)WmL[M1OiMgP+PfQ!Q(R!S."
    "UEV*W1T?U.T/X=LhPfOxMBPQPcQBQrQ|R=RxSBSpS+TxUGXST!X8JGNsJcPG"
    "PpRqQ{R.PtQtMNTEOW00MzOsO?QnQ!QxP)P(Q3R)S]U100UA00U!NPKOPhM?"
    "P-NSPoNbRlPnSp00Rk00TxPgPgPnP,O/PfPiQ)Qw0000S(TxV+00U)MaOoMl"
    "OiM]PaNrQjOtRo00Rw00RqNpNSOV00O(00PlQj0000RkPfQf00Q<000000Sy"
    "00M1MVMfNxOzQiQH0000O/00Of00MfSyMW00Nb00O,000000O(NaP5MpP.00"
    "0000RlJQ00LOLQMR00LP000000Mk00OhM_QhRlMn00NS000000MrOLP-QkRr"
    "000000SqL*000000OdMHP/OWQ%SqMP000000O1O1PXPTQa000000PkMaOVOc"
    "OzP0PoQoNrOcOpP1N9NZNTOQ00NbN*L%")
# Side branches: 2 letters place, 1 letter decay, 1 letter %
BR = (
    "0IdO0Pb}0Qdx0Uq}0ZhB0ad40cd}0ctD0gi}0kh50ld20mdE0nd}0pd.0ptJ"
    "0rd]0rtO0vg-0vik0%d}0(dY0)de0*dt0*tF0+d:0+tR0-d,0-tj0:iA0?d}"
    "0[hK0[dK0]do0_dp0{d}0|dg0|tK0}do0}tG13c}1CdU1DdF1Edo1QdJ1RdM"
    "1SdE1TdV1UdM1Vd}1Xd}1cg}1ciB1ld01mdA1ndK1ntC1oda1otC1pdL1ptH"
    "1xiS1)d01*d51+dU1,dc1,t91-dl1-t61.dW1.tG1/dt1/tL1:tw1:dN1[gc"
    "1]gf1]i425dE26dE27dM29dy2Gg#2Hg92PdA2Qd62RdG2SdY2Tdk2Udd2cge"
    "2dg}2ege2fgh2pdD2qdI2rdP2sdX2!gj2;dK2<dY2=dO2>dc2?dw2[d}2}gA"
    "30gT31gt3EdC3Fdm3GdQ3Hd!3NgA3ScA3adG3bdW3cd+3ddB3kg=3lgi3mgl"
    "3!d13#dk3%d:3(dj3:g33?cJ41d942d}43db44d:45d{45tA46d<46tP4Bg}"
    "4Cgx4Dg}4Tdm4bg34pdO4qdP4yg}4zc}4!g}5KdI5Mdh5Rg}5Sg}5TgF5Ugg"
    "5wg(5xgB5yg25:d15=dB5>dG5{dw61cp61gp62g}64gN65g(6Sda6Zg?6agC"
    "6udO6vdE6wdU6xdQ6ydO6%gc6(g}6)g_6*g}6+gd6,gV76dI78dY7Hg37Pbk"
    "7adC7bdB7cdc7dd:7ed;7kg}7lg}7mg:7ngB7og07-dA7.d97/dG7:d>7<d_"
    "7?gA7[g183e38Dd88EdK8FdU8Gd.8Hds8HtA8Idh8Jd+8JtS8Ng-8kdI8ldO"
    "8mdr8%bg8;d18<d.8=di8>dN8|g}8|q490g}92g49Md39Nd99OdG9QdT9icG"
    "9pdB9qdF9rdM9sdX9tdS9udf9vd=9wd=9xdg9yda9)c.9)qj9*g}9+g}9,g9"
    "9-gB9/g1A7dAA8d9A9dBAAdCABdFACdFADdJATbCAcd9AddIAedHAfdLAgdY"
    "AhdMAidSAjdEAjt1AkdaAldQAsgDA>d1A?d0A[d9A]dAA_dEBPd2BQdABRd9"
    "BSdABTdDBUdGBVdgBWd}Bhg0B*dAB+dFB,dVCGdACHdCCJdECKdNCLdmCXg0"
    "CZgNCvd9CwdACxdKD9d7DAdADBdADPgBDsdDD#g6D%g9D(gAD*gAD<bUD>c3"
    "E6dAE7dEE8dDEBdWEld6Emd7EndAEod9ErdDE#gGE:b9E<cBE>e2F4d1F5d5"
    "F6d9FGgCFKg1FpdCFqdCFrd+FvgAFwg0FygAF/bkF;c4G4d1G5d8G6dAG7dF"
    "G8d}G9d:GEgPGFgYGHg9GmdPGndTGodZGpd)GqdiGugFGup9G,cAG;cAH2dU"
    "H3dQH3t9H4dvH6d}H7dVH7tGHDccHEavHEgAHFgHHFaCHNcXHPcAHgd9HhdB"
    "HidEHqcHHqp9HsaPHsgJHta0Hug8Hui0Hzc}H#cxH#evH-btH/cFH_dFH{dD"
    "H|dIH}dGI0dTI9ciIBa9ICgFIEg2Ifd1IidDIjdFIkdGIoa2IqgHIqi1IsgB"
    "I+bAI-bAI_dAI{dBI|dNI}dMJ0daJ1dXJ2dXJ3dSJ6gSJ6a8J7gNJ8gBJBgX"
    "Jed3JfdAJlcEJ+bgJ?d1J[d9J]dBK3g9KDcJLsbgMUcJMpbaM-gANSg2NbaP"
    "Ngb0NkbPN#gJN,gwN;aiN<aEN=a0ORaUOSaKOebmOyg1OzgFO%cIO(ctO)a4"
    "PBpNPDg1PEg9PHcHPIcsPJa7PSclPYe0PoaIPpcFPqcJPraIPsa4P}cbQ0cj"
    "Q1gNQ4cIQ5cDQ7a8Q8a0QicMQjaqQkahQla6Qma0Q[cbQ]pBQ_cHQ{agR1a1"
    "RIbNRZc0RacLRbcZRcarRdaMReaCRfa1R<aJR=a9R>cER?aeR[coR]aMR_aW"
    "R{a9SGeGSXcCSYcISZcaSacpSbapScaMSdaISeaASfa1Sga3S<pFS=peS=cS"
    "S[aDS]aNS_aAS{aFS|a4S}a7T0aBTGeDTZcATacITbcETccOTdcWTeciTfam"
    "TgaEThaGTia1Tja2T;aJT_amT{amT|aUT}aAU0aBU1a0U2a5U4a2UFbFUccc"
    "UdcrUdg1UeauUfaZUgaMUhaKUia9UjaEU?pZU[ctU_aEU|aDV0aAVKcBVSdA"
    "VTdGVUdgVVdDVccAVecSVgamVhaIViaGVja3Vka3Vwa9V?aIV[cCV|cIV}cv"
    "W0aKW1aCW2a4WJb2WKaiWLaAWdc4WfcFWgcEWhcAWiaqWjcpWkaGWlaJWma9"
    "WnaAWoa0Wpa6WraDWuc4W%bDW)aaW?f7W]cDW_cCW{cIW|cIW}cuX0cbX1aj"
    "X2adX3aCX4aIX5a8X6aHX7a5X8aCX9a1XAaoXIbEXJaGXZcEXacMXccUXdcg"
    "XecaXfaXXgckXhaTXickXjaPXkcCXlaZXvaUX<cGX=cAX>eOX?cDX[cJX]cJ"
    "X_cmX{eSX|apX}c5Y6b3Y7b0YQcMYRcDYTcCYUcFYVcNYWcMYma}Yve9Yyc1"
    "YzcDY#eJY)cAY.c9Y/e9Y:aHY<ePY=a9ZJeIZceAZmcYZneNZoaAZpa4ZqbG"
    "Z.c1Z=eDZ>aSaSamaTeeaUaBaVaAaabKaaa1amaIanaIaoa0apaEazf0a.cH"
    "a;aBa<aCa>a3a}ePb6aSb7fmb7cSb9aQbBaCbDe4bEa9bGe2bLfGbQeSbSaP"
    "bSfCbVf9bVa9bYa4bZa1bba0blfIbofAbpaXbqclbraMbtahb!fBb#a2b*cc"
    "b*fDb,cpb-clb-f9b.aCb/amb:aIb<a1b=a5b>aBb?a4b[eUb{aGc2cIc2fH"
    "c5fFc5e9c6cic7eDc7f0c8afc9eIcAaAcCaKcGf1cJcMcMf0cNaScOcmcPaF"
    "cQaIcSa6cUaFcVaHcVfBcWaNcWf9cbf1cdf2ceffcecAcfcrcgeIcgf1chac"
    "cif4cjeccjfAclcXclfIcpamcpfAcqf9crcacrf0cse2csf0cteNcueNcveB"
    "cwfUcxemc-eRc-f9c/eNc/fGc_fSc{ecc|fEc}cVd1fId1eAd2fQdDfSdDeL"
    "dEfwdFf9dHfcdNfNdPfwdQfQdQeQdRfIdXfGdbfXdkfwdlf9dmfWdpfDdqfw"
    "drfNdtfKdvfSdzfXd-fXd>fwd]fXd}fPe3fKe8fweGfXeNf!eOfPeSfKeUfN"
    "eWfEeXfwedfXeffDekfwemfKenfmeofDepfBeqfEerfwerawesfXeufDewfG"
    "exfSeyfHezfKe#fEe(fS")
# --- DATA END ---

# Where each element starts in MODES (one letter per nucleus)
OFF = [0]
for c in CNT:
    OFF.append(OFF[-1] + c)
SYM = SYM.split()

# letter: (name in table, change of Z, change of N)
MOVE = {
    "a": ("alpha", -2, -2),
    "b": ("beta-", 1, -1),
    "c": ("beta+", -1, 1),
    "e": ("EC", -1, 1),
    "p": ("p", -1, 0),
    "q": ("2p", -2, 0),
    "n": ("n", 0, -1),
    "m": ("2n", 0, -2),
    "d": ("b-,n", 1, -2),
    "t": ("b-,2n", 1, -3),
    "g": ("b+,p", -2, 1),
    "h": ("b-,a", -1, -3),
    "i": ("b+,a", -3, -1),
}
WORD = {"alpha": "alpha", "beta-": "beta minus", "beta+": "beta plus",
        "EC": "electron capture", "p": "proton out", "2p": "2 protons out",
        "n": "neutron out", "2n": "2 neutrons out", "b-,n": "beta- + neutron",
        "b-,2n": "beta- + 2 neutrons", "b+,p": "beta+ + proton",
        "b-,a": "beta- + alpha", "b+,a": "beta+ + alpha"}
# Half-life units: (name, seconds, smallest k, biggest k)
# A half-life is m * 10^(k-1) units, m = 10..99 (2 digits)
UNITS = [("s", 1, -24, -10), ("ns", 1, 0, 2), ("us", 1, 0, 2),
         ("ms", 1, 0, 2), ("s", 1, 0, 1), ("min", 1, 0, 1),
         ("h", 1, 0, 1), ("d", 1, 0, 2), ("y", 1, 0, 32)]
ABC = ("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
       "abcdefghijklmnopqrstuvwxyz!#%()*+,-./:;<=>?[]_{|}")
LINES = 9          # lines per screen

def find(z, n):    # place in MODES, or -1 if the nucleus is not known
    if z < 0 or z >= len(CNT) or n < NMIN[z] or n >= NMIN[z] + CNT[z]:
        return -1
    return OFF[z] + n - NMIN[z]

def mode(z, n):
    i = find(z, n)
    return "x" if i < 0 else MODES[i]

def name(z, n):
    return SYM[z] + "-" + str(z + n)

def pad(s, w):     # MicroPython has no ljust
    return s + " " * (w - len(s))

def half(z, n):    # half-life as short text, "" if not known
    i = find(z, n)
    if i < 0:
        return ""
    v = ABC.find(HL[2 * i]) * 85 + ABC.find(HL[2 * i + 1])
    if v == 0:
        return ""
    v -= 1
    for u in UNITS:
        size = (u[3] - u[2] + 1) * 90
        if v < size:
            k = u[2] + v // 90
            m = v % 90 + 10          # 2 digits: 10..99
            if k < 0 or k > 5:
                t = str(m // 10) + ("." + str(m % 10) if m % 10 else "")
                return t + "E" + str(k) + " " + u[0]
            if k == 0:
                t = str(m // 10) + ("." + str(m % 10) if m % 10 else "")
            else:
                t = str(m * 10 ** (k - 1))
            return t + " " + u[0]
        v -= size
    return ""

def num2(s):       # 2 letters of ABC -> number
    return ABC.find(s[0]) * 85 + ABC.find(s[1])

def sides(z, n):   # other ways to decay: list of (letter, %)
    i = find(z, n)
    out = []
    if i < 0:
        return out
    lo, hi = 0, len(BR) // 4
    while lo < hi:               # binary search, BR is sorted by place
        mid = (lo + hi) // 2
        if num2(BR[4 * mid:4 * mid + 2]) < i:
            lo = mid + 1
        else:
            hi = mid
    while 4 * lo < len(BR) and num2(BR[4 * lo:4 * lo + 2]) == i:
        c = ABC.find(BR[4 * lo + 3])
        out.append((BR[4 * lo + 2], (c + 1) / 10 if c < 10 else c - 8))
        lo += 1
    return out

def pct(p):        # percent as short text
    if p >= 99.95:
        return "100%"
    if p >= 9.95:
        return str(int(p + 0.5)) + "%"
    if p >= 0.995:
        t = str(int(p * 10 + 0.5) / 10)
        return (t[:-2] if t[-2:] == ".0" else t) + "%"
    if p >= 0.01:
        return str(int(p * 100 + 0.5) / 100) + "%"
    return "<0.01%"

def ends(z, n):    # all possible end nuclei with their chance in %
    todo = {(z, n): 100.0}
    res = {}
    loops = 0
    while todo and loops < 3000:
        loops += 1
        best = None              # heaviest first: all paths meet there
        for k in todo:
            if best is None or k[0] + k[1] > best[0] + best[1]:
                best = k
        p = todo.pop(best)
        m = mode(best[0], best[1])
        if m not in MOVE:
            res[best] = res.get(best, 0) + p
            continue
        go = [(m, 100.0)]
        for s in sides(best[0], best[1]):
            go.append(s)
            go[0] = (m, go[0][1] - s[1])
        for g in go:
            if g[0] not in MOVE:         # side branch = fission
                k = best + ("f",)
                res[k] = res.get(k, 0) + p * g[1] / 100
                continue
            w = MOVE[g[0]]
            k = (best[0] + w[1], best[1] + w[2])
            todo[k] = todo.get(k, 0) + p * max(g[1], 0) / 100
    out = [(res[k], k) for k in res if res[k] >= 0.005]
    out.sort()
    out.reverse()
    return out

class Screen:      # prints lines and waits when the screen is full
    def __init__(self):
        self.used = 0
    def put(self, s):
        if self.used >= LINES:
            self.wait()
        print(s)
        self.used += 1
    def wait(self):
        print("[any key]")
        if wait_key:
            wait_key()
            disp_clr()
        else:
            input("")          # on the PC: press enter
        self.used = 0

def ask(text):
    try:
        return input(text).strip()
    except:
        return ""

def read_z(s):     # a number or an element symbol (U, Pb, ...)
    for i in range(len(SYM)):
        if SYM[i].lower() == s.lower():
            return i
    return int(s)

def solve(z, n):
    sc = Screen()
    if disp_clr:
        disp_clr()
    if mode(z, n) == "x":
        sc.put("Z=" + str(z) + " N=" + str(n))
        sc.put("Not on the nuclide chart.")
        sc.put("No such nucleus is known.")
        return
    z0, n0 = z, n
    sc.put(name(z, n) + "  Z=" + str(z) + " N=" + str(n))
    sc.put("Nucleus Decay T1/2")
    count = {}
    steps = 0
    while steps < 80:
        m = mode(z, n)
        if m in MOVE:
            w = MOVE[m]
            sc.put(pad(name(z, n), 7) + pad(w[0], 6) + half(z, n))
            for s in sides(z, n):
                t = MOVE[s[0]][0] if s[0] in MOVE else "fission"
                sc.put("   or " + t + " " + pct(s[1]))
            count[w[0]] = count.get(w[0], 0) + 1
            z, n = z + w[1], n + w[2]
            steps += 1
            continue
        if m == "s":
            sc.put(pad(name(z, n), 7) + "stable")
        elif m == "l":
            sc.put(pad(name(z, n), 7) + "~stable")
        elif m == "f":
            sc.put(pad(name(z, n), 7) + "fission")
        else:
            sc.put(pad(name(z, n), 7) + "?")
        break
    sc.wait()
    sc.put("Result: " + name(z, n))
    sc.put("Z=" + str(z) + " N=" + str(n) + " A=" + str(z + n))
    m = mode(z, n)
    if m == "s":
        sc.put("Stable.")
    elif m == "l":
        sc.put("Practically stable:")
        sc.put("T1/2 = " + half(z, n))
    elif m == "f":
        sc.put("Splits in two (fission).")
        sc.put("The pieces vary, so")
        sc.put("there is no single end.")
    elif m == "x":
        sc.put("Not on the nuclide chart.")
    else:
        sc.put("Decay not known.")
    if steps == 0:
        sc.put("No decay steps.")
        return
    sc.put("Main path: " + str(steps) + (" step" if steps == 1 else " steps"))
    for k in count:
        sc.put(str(count[k]) + "x " + WORD[k])
    if gc:
        gc.collect()
    e = ends(z0, n0)
    if len(e) < 2:
        return
    sc.wait()
    sc.put("All possible ends:")
    for p, k in e[:12]:
        t = pad(name(k[0], k[1]), 7) + pad(pct(p), 7)
        m = "f" if len(k) > 2 else mode(k[0], k[1])
        if m == "l":
            t += "~stable"
        elif m == "f":
            t += "fission"
        elif m != "s":
            t += "?"
        sc.put(t)
    if len(e) > 12:
        sc.put("+ " + str(len(e) - 12) + " more (small)")

print("Decay chain solver")
print("Type protons and neutrons.")
print("Empty = quit.")
while True:
    s = ask("Protons Z? ")
    if not s:
        break
    try:
        z = read_z(s)
        n = int(ask("Neutrons N? "))
    except:
        print("Please type a number.")
        continue
    solve(z, n)
    if gc:
        gc.collect()
    print("")
