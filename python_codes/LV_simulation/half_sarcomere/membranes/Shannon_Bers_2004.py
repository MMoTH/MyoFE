# Size of variable arrays:
sizeAlgebraic = 111
sizeStates = 39
sizeConstants = 137
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (millisecond)"
    legend_states[0] = "V in component cell (millivolt)"
    legend_algebraic[36] = "i_Na in component INa (microA_per_microF)"
    legend_algebraic[39] = "i_Nab in component INab (microA_per_microF)"
    legend_algebraic[43] = "i_NaK in component INaK (microA_per_microF)"
    legend_algebraic[45] = "i_Kr in component IKr (microA_per_microF)"
    legend_algebraic[53] = "i_Ks in component IKs (microA_per_microF)"
    legend_algebraic[54] = "i_Kp in component IKp (microA_per_microF)"
    legend_algebraic[55] = "i_tos in component Itos (microA_per_microF)"
    legend_algebraic[56] = "i_tof in component Itof (microA_per_microF)"
    legend_algebraic[60] = "i_K1 in component IK1 (microA_per_microF)"
    legend_algebraic[79] = "i_NaCa in component INaCa (microA_per_microF)"
    legend_algebraic[61] = "i_Cl_Ca in component ICl_Ca (microA_per_microF)"
    legend_algebraic[62] = "i_Clb in component IClb (microA_per_microF)"
    legend_algebraic[71] = "i_CaL in component ICaL (microA_per_microF)"
    legend_algebraic[88] = "i_Cab in component ICab (microA_per_microF)"
    legend_algebraic[85] = "i_Cap in component ICap (microA_per_microF)"
    legend_algebraic[4] = "i_Stim in component cell (microA_per_microF)"
    legend_constants[0] = "stim_start in component cell (millisecond)"
    legend_constants[1] = "stim_period in component cell (millisecond)"
    legend_constants[2] = "stim_duration in component cell (millisecond)"
    legend_constants[3] = "stim_amplitude in component cell (microA_per_microF)"
    legend_constants[4] = "Ko in component model_parameters (millimolar)"
    legend_constants[5] = "Nao in component model_parameters (millimolar)"
    legend_constants[6] = "Cao in component model_parameters (millimolar)"
    legend_constants[7] = "Clo in component model_parameters (millimolar)"
    legend_constants[8] = "Ki in component model_parameters (millimolar)"
    legend_constants[9] = "Mgi in component model_parameters (millimolar)"
    legend_constants[10] = "Cli in component model_parameters (millimolar)"
    legend_constants[11] = "R in component model_parameters (joule_per_kilomole_kelvin)"
    legend_constants[12] = "T in component model_parameters (kelvin)"
    legend_constants[13] = "F in component model_parameters (coulomb_per_mole)"
    legend_constants[14] = "Cm in component model_parameters (farad)"
    legend_constants[15] = "cell_length in component model_parameters (micrometre)"
    legend_constants[16] = "cell_radius in component model_parameters (micrometre)"
    legend_constants[132] = "Vol_Cell in component model_parameters (litre)"
    legend_constants[133] = "Vol_SR in component model_parameters (litre)"
    legend_constants[134] = "Vol_SL in component model_parameters (litre)"
    legend_constants[135] = "Vol_jct in component model_parameters (litre)"
    legend_constants[136] = "Vol_myo in component model_parameters (litre)"
    legend_algebraic[18] = "E_Na_jct in component reversal_potentials (millivolt)"
    legend_algebraic[28] = "E_Na_SL in component reversal_potentials (millivolt)"
    legend_algebraic[30] = "E_Ca_jct in component reversal_potentials (millivolt)"
    legend_algebraic[32] = "E_Ca_SL in component reversal_potentials (millivolt)"
    legend_constants[123] = "E_K in component reversal_potentials (millivolt)"
    legend_constants[124] = "E_Cl in component reversal_potentials (millivolt)"
    legend_states[1] = "Nai in component Na_buffer (millimolar)"
    legend_states[2] = "Na_jct in component Na_buffer (millimolar)"
    legend_states[3] = "Na_SL in component Na_buffer (millimolar)"
    legend_states[4] = "Cai in component Ca_buffer (millimolar)"
    legend_states[5] = "Ca_jct in component Ca_buffer (millimolar)"
    legend_states[6] = "Ca_SL in component Ca_buffer (millimolar)"
    legend_algebraic[34] = "i_Na_jct in component INa (microA_per_microF)"
    legend_algebraic[35] = "i_Na_SL in component INa (microA_per_microF)"
    legend_constants[17] = "G_INa in component INa (milliS_per_microF)"
    legend_constants[18] = "Fx_Na_jct in component INa (dimensionless)"
    legend_constants[19] = "Fx_Na_SL in component INa (dimensionless)"
    legend_states[7] = "m in component INa_m_gate (dimensionless)"
    legend_states[8] = "h in component INa_h_gate (dimensionless)"
    legend_states[9] = "j in component INa_j_gate (dimensionless)"
    legend_algebraic[33] = "openProb in component INa (dimensionless)"
    legend_algebraic[0] = "alpha_h in component INa_h_gate (per_millisecond)"
    legend_algebraic[14] = "beta_h in component INa_h_gate (per_millisecond)"
    legend_algebraic[1] = "alpha_j in component INa_j_gate (per_millisecond)"
    legend_algebraic[15] = "beta_j in component INa_j_gate (per_millisecond)"
    legend_algebraic[2] = "alpha_m in component INa_m_gate (per_millisecond)"
    legend_algebraic[16] = "beta_m in component INa_m_gate (per_millisecond)"
    legend_algebraic[37] = "i_Nab_jct in component INab (microA_per_microF)"
    legend_algebraic[38] = "i_Nab_SL in component INab (microA_per_microF)"
    legend_constants[20] = "G_NaBk in component INab (milliS_per_microF)"
    legend_constants[21] = "Fx_NaBk_jct in component INab (dimensionless)"
    legend_constants[22] = "Fx_NaBk_SL in component INab (dimensionless)"
    legend_algebraic[41] = "i_NaK_jct in component INaK (microA_per_microF)"
    legend_algebraic[42] = "i_NaK_SL in component INaK (microA_per_microF)"
    legend_algebraic[40] = "f_NaK in component INaK (dimensionless)"
    legend_constants[23] = "H_NaK in component INaK (dimensionless)"
    legend_constants[24] = "Km_Nai in component INaK (millimolar)"
    legend_constants[25] = "Km_Ko in component INaK (millimolar)"
    legend_constants[26] = "I_NaK_max in component INaK (microA_per_microF)"
    legend_constants[125] = "sigma in component INaK (dimensionless)"
    legend_constants[27] = "Fx_NaK_jct in component INaK (dimensionless)"
    legend_constants[28] = "Fx_NaK_SL in component INaK (dimensionless)"
    legend_constants[126] = "G_IKr in component IKr (milliS_per_microF)"
    legend_states[10] = "Xr in component IKr_Xr_gate (dimensionless)"
    legend_algebraic[44] = "Rr in component IKr_Rr_gate (dimensionless)"
    legend_algebraic[3] = "Xr_infinity in component IKr_Xr_gate (dimensionless)"
    legend_algebraic[17] = "tau_Xr in component IKr_Xr_gate (millisecond)"
    legend_algebraic[51] = "i_Ks_jct in component IKs (microA_per_microF)"
    legend_algebraic[52] = "i_Ks_SL in component IKs (microA_per_microF)"
    legend_algebraic[49] = "G_Ks_SL in component IKs (milliS_per_microF)"
    legend_algebraic[48] = "G_Ks_jct in component IKs (milliS_per_microF)"
    legend_constants[29] = "Fx_Ks_jct in component IKs (dimensionless)"
    legend_constants[30] = "Fx_Ks_SL in component IKs (dimensionless)"
    legend_algebraic[50] = "E_Ks in component IKs (millivolt)"
    legend_constants[31] = "pKNa in component IKs (dimensionless)"
    legend_algebraic[46] = "pCa_jct in component IKs (dimensionless)"
    legend_algebraic[47] = "pCa_SL in component IKs (dimensionless)"
    legend_states[11] = "Xs in component IKs_Xs_gate (dimensionless)"
    legend_algebraic[5] = "Xs_infinity in component IKs_Xs_gate (dimensionless)"
    legend_algebraic[19] = "tau_Xs in component IKs_Xs_gate (millisecond)"
    legend_constants[32] = "g_Kp in component IKp (milliS_per_microF)"
    legend_constants[33] = "G_tos in component Itos (milliS_per_microF)"
    legend_states[12] = "Y_tos in component Itos_Y_gate (dimensionless)"
    legend_states[13] = "X_tos in component Itos_X_gate (dimensionless)"
    legend_states[14] = "R_tos in component Itos_R_gate (dimensionless)"
    legend_algebraic[6] = "X_tos_infinity in component Itos_X_gate (dimensionless)"
    legend_algebraic[20] = "tau_X_tos in component Itos_X_gate (millisecond)"
    legend_algebraic[7] = "Y_tos_infinity in component Itos_Y_gate (dimensionless)"
    legend_algebraic[21] = "tau_Y_tos in component Itos_Y_gate (millisecond)"
    legend_algebraic[8] = "R_tos_infinity in component Itos_R_gate (dimensionless)"
    legend_algebraic[22] = "tau_R_tos in component Itos_R_gate (millisecond)"
    legend_constants[34] = "G_tof in component Itof (milliS_per_microF)"
    legend_states[15] = "Y_tof in component Itof_Y_gate (dimensionless)"
    legend_states[16] = "X_tof in component Itof_X_gate (dimensionless)"
    legend_algebraic[9] = "X_tof_infinity in component Itof_X_gate (dimensionless)"
    legend_algebraic[23] = "tau_X_tof in component Itof_X_gate (millisecond)"
    legend_algebraic[10] = "Y_tof_infinity in component Itof_Y_gate (dimensionless)"
    legend_algebraic[24] = "tau_Y_tof in component Itof_Y_gate (millisecond)"
    legend_constants[127] = "G_K1 in component IK1 (milliS_per_microF)"
    legend_algebraic[59] = "K1_infinity in component IK1_K1_gate (dimensionless)"
    legend_algebraic[57] = "alpha_K1 in component IK1_K1_gate (per_millisecond)"
    legend_algebraic[58] = "beta_K1 in component IK1_K1_gate (per_millisecond)"
    legend_constants[35] = "G_Cl in component ICl_Ca (milliS_per_microF)"
    legend_constants[36] = "Kd_ClCa in component ICl_Ca (millimolar)"
    legend_constants[37] = "Fx_Cl_jct in component ICl_Ca (dimensionless)"
    legend_constants[38] = "Fx_Cl_SL in component ICl_Ca (dimensionless)"
    legend_constants[39] = "G_ClBk in component IClb (milliS_per_microF)"
    legend_algebraic[65] = "i_CaL_Ca_SL in component ICaL (microA_per_microF)"
    legend_algebraic[68] = "i_CaL_Ca_jct in component ICaL (microA_per_microF)"
    legend_algebraic[66] = "i_CaL_Na_SL in component ICaL (microA_per_microF)"
    legend_algebraic[69] = "i_CaL_Na_jct in component ICaL (microA_per_microF)"
    legend_algebraic[70] = "i_CaL_K in component ICaL (microA_per_microF)"
    legend_constants[40] = "PCa in component ICaL (litre_per_farad_millisecond)"
    legend_constants[41] = "PNa in component ICaL (litre_per_farad_millisecond)"
    legend_constants[42] = "PK in component ICaL (litre_per_farad_millisecond)"
    legend_constants[43] = "Fx_ICaL_jct in component ICaL (dimensionless)"
    legend_constants[44] = "Fx_ICaL_SL in component ICaL (dimensionless)"
    legend_constants[45] = "gamma_Cai in component ICaL (dimensionless)"
    legend_constants[46] = "gamma_Cao in component ICaL (dimensionless)"
    legend_constants[47] = "gamma_Nai in component ICaL (dimensionless)"
    legend_constants[48] = "gamma_Nao in component ICaL (dimensionless)"
    legend_constants[49] = "gamma_Ki in component ICaL (dimensionless)"
    legend_constants[50] = "gamma_Ko in component ICaL (dimensionless)"
    legend_constants[51] = "Q10_CaL in component ICaL (dimensionless)"
    legend_constants[128] = "Q_CaL in component ICaL (dimensionless)"
    legend_states[17] = "d in component ICaL_d_gate (dimensionless)"
    legend_states[18] = "f in component ICaL_f_gate (dimensionless)"
    legend_algebraic[64] = "fCa_SL in component ICaL_fCa_gate (dimensionless)"
    legend_algebraic[67] = "fCa_jct in component ICaL_fCa_gate (dimensionless)"
    legend_algebraic[63] = "temp in component ICaL (coulomb_per_mole)"
    legend_algebraic[11] = "d_infinity in component ICaL_d_gate (dimensionless)"
    legend_algebraic[25] = "tau_d in component ICaL_d_gate (millisecond)"
    legend_algebraic[12] = "f_infinity in component ICaL_f_gate (dimensionless)"
    legend_algebraic[26] = "tau_f in component ICaL_f_gate (millisecond)"
    legend_states[19] = "fCaB_SL in component ICaL_fCa_gate (dimensionless)"
    legend_states[20] = "fCaB_jct in component ICaL_fCa_gate (dimensionless)"
    legend_algebraic[77] = "i_NaCa_SL in component INaCa (microA_per_microF)"
    legend_algebraic[76] = "i_NaCa_jct in component INaCa (microA_per_microF)"
    legend_constants[52] = "V_max in component INaCa (microA_per_microF)"
    legend_constants[53] = "Fx_NCX_jct in component INaCa (dimensionless)"
    legend_constants[54] = "Fx_NCX_SL in component INaCa (dimensionless)"
    legend_constants[55] = "Q10_NCX in component INaCa (dimensionless)"
    legend_constants[129] = "Q_NCX in component INaCa (dimensionless)"
    legend_constants[56] = "K_mNai in component INaCa (millimolar)"
    legend_constants[57] = "K_mCao in component INaCa (millimolar)"
    legend_constants[58] = "K_mNao in component INaCa (millimolar)"
    legend_constants[59] = "K_mCai in component INaCa (millimolar)"
    legend_algebraic[74] = "Ka_SL in component INaCa (dimensionless)"
    legend_algebraic[75] = "Ka_jct in component INaCa (dimensionless)"
    legend_constants[60] = "Kd_act in component INaCa (millimolar)"
    legend_constants[61] = "ksat in component INaCa (dimensionless)"
    legend_constants[62] = "eta in component INaCa (dimensionless)"
    legend_constants[63] = "HNa in component INaCa (dimensionless)"
    legend_algebraic[72] = "temp_jct in component INaCa (millimolar4)"
    legend_algebraic[73] = "temp_SL in component INaCa (millimolar4)"
    legend_algebraic[84] = "i_Cap_SL in component ICap (microA_per_microF)"
    legend_algebraic[82] = "i_Cap_jct in component ICap (microA_per_microF)"
    legend_constants[64] = "Fx_SLCaP_jct in component ICap (dimensionless)"
    legend_constants[65] = "Fx_SLCaP_SL in component ICap (dimensionless)"
    legend_constants[66] = "Q10_SLCaP in component ICap (dimensionless)"
    legend_constants[130] = "Q_SLCaP in component ICap (dimensionless)"
    legend_constants[67] = "Km in component ICap (millimolar)"
    legend_constants[68] = "H in component ICap (dimensionless)"
    legend_constants[69] = "V_maxAF in component ICap (microA_per_microF)"
    legend_algebraic[86] = "i_Cab_jct in component ICab (microA_per_microF)"
    legend_algebraic[87] = "i_Cab_SL in component ICab (microA_per_microF)"
    legend_constants[70] = "G_CaBk in component ICab (milliS_per_microF)"
    legend_constants[71] = "Fx_CaBk_jct in component ICab (dimensionless)"
    legend_constants[72] = "Fx_CaBk_SL in component ICab (dimensionless)"
    legend_algebraic[89] = "j_rel_SR in component Jrel_SR (millimolar_per_millisecond)"
    legend_constants[73] = "Max_SR in component Jrel_SR (dimensionless)"
    legend_constants[74] = "Min_SR in component Jrel_SR (dimensionless)"
    legend_constants[75] = "EC50_SR in component Jrel_SR (millimolar)"
    legend_states[21] = "R in component Jrel_SR (dimensionless)"
    legend_states[22] = "I in component Jrel_SR (dimensionless)"
    legend_algebraic[31] = "RI in component Jrel_SR (dimensionless)"
    legend_states[23] = "O in component Jrel_SR (dimensionless)"
    legend_constants[76] = "ks in component Jrel_SR (per_millisecond)"
    legend_constants[77] = "koCa in component Jrel_SR (per_millimolar2_per_millisecond)"
    legend_constants[78] = "kom in component Jrel_SR (per_millisecond)"
    legend_constants[79] = "kiCa in component Jrel_SR (per_millimolar_per_millisecond)"
    legend_constants[80] = "kim in component Jrel_SR (per_millisecond)"
    legend_constants[81] = "HSR in component Jrel_SR (dimensionless)"
    legend_states[24] = "Ca_SR in component Ca_buffer (millimolar)"
    legend_algebraic[13] = "kCaSR in component Jrel_SR (dimensionless)"
    legend_algebraic[27] = "koSRCa in component Jrel_SR (per_millimolar2_per_millisecond)"
    legend_algebraic[29] = "kiSRCa in component Jrel_SR (per_millimolar_per_millisecond)"
    legend_algebraic[90] = "j_leak_SR in component Jleak_SR (millimolar_per_millisecond)"
    legend_constants[82] = "KSRleak in component Jleak_SR (per_millisecond)"
    legend_algebraic[91] = "j_pump_SR in component Jpump_SR (millimolar_per_millisecond)"
    legend_constants[83] = "V_max in component Jpump_SR (millimolar_per_millisecond)"
    legend_constants[84] = "Q10_SRCaP in component Jpump_SR (dimensionless)"
    legend_constants[131] = "Q_SRCaP in component Jpump_SR (dimensionless)"
    legend_constants[85] = "Kmf in component Jpump_SR (millimolar)"
    legend_constants[86] = "Kmr in component Jpump_SR (millimolar)"
    legend_constants[87] = "H in component Jpump_SR (dimensionless)"
    legend_algebraic[78] = "J_Na_jct_SL in component ion_diffusion (millimole_per_millisecond)"
    legend_algebraic[81] = "J_Na_SL_myo in component ion_diffusion (millimole_per_millisecond)"
    legend_algebraic[93] = "J_Ca_jct_SL in component ion_diffusion (millimole_per_millisecond)"
    legend_algebraic[94] = "J_Ca_SL_myo in component ion_diffusion (millimole_per_millisecond)"
    legend_states[25] = "Na_SL_buf in component Na_buffer (millimolar)"
    legend_states[26] = "Na_jct_buf in component Na_buffer (millimolar)"
    legend_constants[88] = "Bmax_SL in component Na_buffer (millimolar)"
    legend_constants[89] = "Bmax_jct in component Na_buffer (millimolar)"
    legend_constants[90] = "kon in component Na_buffer (per_millimolar_per_millisecond)"
    legend_constants[91] = "koff in component Na_buffer (per_millisecond)"
    legend_algebraic[80] = "dNa_jct_buf in component Na_buffer (millimolar_per_millisecond)"
    legend_algebraic[83] = "dNa_SL_buf in component Na_buffer (millimolar_per_millisecond)"
    legend_states[27] = "Ca_SLB_SL in component Ca_buffer (millimolar)"
    legend_states[28] = "Ca_SLB_jct in component Ca_buffer (millimolar)"
    legend_states[29] = "Ca_SLHigh_SL in component Ca_buffer (millimolar)"
    legend_states[30] = "Ca_SLHigh_jct in component Ca_buffer (millimolar)"
    legend_states[31] = "Ca_Calsequestrin in component Ca_buffer (millimolar)"
    legend_constants[92] = "Bmax_SLB_SL in component Ca_buffer (millimolar)"
    legend_constants[93] = "Bmax_SLB_jct in component Ca_buffer (millimolar)"
    legend_constants[94] = "Bmax_SLHigh_SL in component Ca_buffer (millimolar)"
    legend_constants[95] = "Bmax_SLHigh_jct in component Ca_buffer (millimolar)"
    legend_constants[96] = "Bmax_Calsequestrin in component Ca_buffer (millimolar)"
    legend_constants[97] = "kon_SL in component Ca_buffer (per_millimolar_per_millisecond)"
    legend_constants[98] = "kon_Calsequestrin in component Ca_buffer (per_millimolar_per_millisecond)"
    legend_constants[99] = "koff_SLB in component Ca_buffer (per_millisecond)"
    legend_constants[100] = "koff_SLHigh in component Ca_buffer (per_millisecond)"
    legend_constants[101] = "koff_Calsequestrin in component Ca_buffer (per_millisecond)"
    legend_algebraic[92] = "dCalsequestrin in component Ca_buffer (millimolar_per_millisecond)"
    legend_algebraic[110] = "dCa_cytosol_tot_bound in component cytosolic_Ca_buffer (millimolar_per_millisecond)"
    legend_algebraic[97] = "dCa_SLB_SL in component Ca_buffer (millimolar_per_millisecond)"
    legend_algebraic[95] = "dCa_SLB_jct in component Ca_buffer (millimolar_per_millisecond)"
    legend_algebraic[100] = "dCa_SLHigh_SL in component Ca_buffer (millimolar_per_millisecond)"
    legend_algebraic[98] = "dCa_SLHigh_jct in component Ca_buffer (millimolar_per_millisecond)"
    legend_algebraic[101] = "dCa_jct_tot_bound in component Ca_buffer (millimolar_per_millisecond)"
    legend_algebraic[104] = "dCa_SL_tot_bound in component Ca_buffer (millimolar_per_millisecond)"
    legend_algebraic[103] = "i_Ca_jct_tot in component Ca_buffer (microA_per_microF)"
    legend_algebraic[106] = "i_Ca_SL_tot in component Ca_buffer (microA_per_microF)"
    legend_states[32] = "Ca_TroponinC in component cytosolic_Ca_buffer (millimolar)"
    legend_states[33] = "Ca_TroponinC_Ca_Mg in component cytosolic_Ca_buffer (millimolar)"
    legend_states[34] = "Mg_TroponinC_Ca_Mg in component cytosolic_Ca_buffer (millimolar)"
    legend_states[35] = "Ca_Calmodulin in component cytosolic_Ca_buffer (millimolar)"
    legend_states[36] = "Ca_Myosin in component cytosolic_Ca_buffer (millimolar)"
    legend_states[37] = "Mg_Myosin in component cytosolic_Ca_buffer (millimolar)"
    legend_states[38] = "Ca_SRB in component cytosolic_Ca_buffer (millimolar)"
    legend_constants[102] = "Bmax_TroponinC in component cytosolic_Ca_buffer (millimolar)"
    legend_constants[103] = "Bmax_TroponinC_Ca_Mg_Ca in component cytosolic_Ca_buffer (millimolar)"
    legend_constants[104] = "Bmax_TroponinC_Ca_Mg_Mg in component cytosolic_Ca_buffer (millimolar)"
    legend_constants[105] = "Bmax_Calmodulin in component cytosolic_Ca_buffer (millimolar)"
    legend_constants[106] = "Bmax_Myosin_Ca in component cytosolic_Ca_buffer (millimolar)"
    legend_constants[107] = "Bmax_Myosin_Mg in component cytosolic_Ca_buffer (millimolar)"
    legend_constants[108] = "Bmax_SRB in component cytosolic_Ca_buffer (millimolar)"
    legend_constants[109] = "kon_TroponinC in component cytosolic_Ca_buffer (per_millimolar_per_millisecond)"
    legend_constants[110] = "kon_TroponinC_Ca_Mg_Ca in component cytosolic_Ca_buffer (per_millimolar_per_millisecond)"
    legend_constants[111] = "kon_TroponinC_Ca_Mg_Mg in component cytosolic_Ca_buffer (per_millimolar_per_millisecond)"
    legend_constants[112] = "kon_Calmodulin in component cytosolic_Ca_buffer (per_millimolar_per_millisecond)"
    legend_constants[113] = "kon_Myosin_Ca in component cytosolic_Ca_buffer (per_millimolar_per_millisecond)"
    legend_constants[114] = "kon_Myosin_Mg in component cytosolic_Ca_buffer (per_millimolar_per_millisecond)"
    legend_constants[115] = "kon_SRB in component cytosolic_Ca_buffer (per_millimolar_per_millisecond)"
    legend_constants[116] = "koff_TroponinC in component cytosolic_Ca_buffer (per_millisecond)"
    legend_constants[117] = "koff_TroponinC_Ca_Mg_Ca in component cytosolic_Ca_buffer (per_millisecond)"
    legend_constants[118] = "koff_TroponinC_Ca_Mg_Mg in component cytosolic_Ca_buffer (per_millisecond)"
    legend_constants[119] = "koff_Calmodulin in component cytosolic_Ca_buffer (per_millisecond)"
    legend_constants[120] = "koff_Myosin_Ca in component cytosolic_Ca_buffer (per_millisecond)"
    legend_constants[121] = "koff_Myosin_Mg in component cytosolic_Ca_buffer (per_millisecond)"
    legend_constants[122] = "koff_SRB in component cytosolic_Ca_buffer (per_millisecond)"
    legend_algebraic[96] = "dCa_TroponinC in component cytosolic_Ca_buffer (millimolar_per_millisecond)"
    legend_algebraic[99] = "dCa_TroponinC_Ca_Mg in component cytosolic_Ca_buffer (millimolar_per_millisecond)"
    legend_algebraic[102] = "dMg_TroponinC_Ca_Mg in component cytosolic_Ca_buffer (millimolar_per_millisecond)"
    legend_algebraic[105] = "dCa_Calmodulin in component cytosolic_Ca_buffer (millimolar_per_millisecond)"
    legend_algebraic[107] = "dCa_Myosin in component cytosolic_Ca_buffer (millimolar_per_millisecond)"
    legend_algebraic[108] = "dMg_Myosin in component cytosolic_Ca_buffer (millimolar_per_millisecond)"
    legend_algebraic[109] = "dCa_SRB in component cytosolic_Ca_buffer (millimolar_per_millisecond)"
    legend_rates[0] = "d/dt V in component cell (millivolt)"
    legend_rates[8] = "d/dt h in component INa_h_gate (dimensionless)"
    legend_rates[9] = "d/dt j in component INa_j_gate (dimensionless)"
    legend_rates[7] = "d/dt m in component INa_m_gate (dimensionless)"
    legend_rates[10] = "d/dt Xr in component IKr_Xr_gate (dimensionless)"
    legend_rates[11] = "d/dt Xs in component IKs_Xs_gate (dimensionless)"
    legend_rates[13] = "d/dt X_tos in component Itos_X_gate (dimensionless)"
    legend_rates[12] = "d/dt Y_tos in component Itos_Y_gate (dimensionless)"
    legend_rates[14] = "d/dt R_tos in component Itos_R_gate (dimensionless)"
    legend_rates[16] = "d/dt X_tof in component Itof_X_gate (dimensionless)"
    legend_rates[15] = "d/dt Y_tof in component Itof_Y_gate (dimensionless)"
    legend_rates[17] = "d/dt d in component ICaL_d_gate (dimensionless)"
    legend_rates[18] = "d/dt f in component ICaL_f_gate (dimensionless)"
    legend_rates[19] = "d/dt fCaB_SL in component ICaL_fCa_gate (dimensionless)"
    legend_rates[20] = "d/dt fCaB_jct in component ICaL_fCa_gate (dimensionless)"
    legend_rates[21] = "d/dt R in component Jrel_SR (dimensionless)"
    legend_rates[23] = "d/dt O in component Jrel_SR (dimensionless)"
    legend_rates[22] = "d/dt I in component Jrel_SR (dimensionless)"
    legend_rates[26] = "d/dt Na_jct_buf in component Na_buffer (millimolar)"
    legend_rates[25] = "d/dt Na_SL_buf in component Na_buffer (millimolar)"
    legend_rates[2] = "d/dt Na_jct in component Na_buffer (millimolar)"
    legend_rates[3] = "d/dt Na_SL in component Na_buffer (millimolar)"
    legend_rates[1] = "d/dt Nai in component Na_buffer (millimolar)"
    legend_rates[31] = "d/dt Ca_Calsequestrin in component Ca_buffer (millimolar)"
    legend_rates[27] = "d/dt Ca_SLB_SL in component Ca_buffer (millimolar)"
    legend_rates[28] = "d/dt Ca_SLB_jct in component Ca_buffer (millimolar)"
    legend_rates[29] = "d/dt Ca_SLHigh_SL in component Ca_buffer (millimolar)"
    legend_rates[30] = "d/dt Ca_SLHigh_jct in component Ca_buffer (millimolar)"
    legend_rates[24] = "d/dt Ca_SR in component Ca_buffer (millimolar)"
    legend_rates[5] = "d/dt Ca_jct in component Ca_buffer (millimolar)"
    legend_rates[6] = "d/dt Ca_SL in component Ca_buffer (millimolar)"
    legend_rates[4] = "d/dt Cai in component Ca_buffer (millimolar)"
    legend_rates[32] = "d/dt Ca_TroponinC in component cytosolic_Ca_buffer (millimolar)"
    legend_rates[33] = "d/dt Ca_TroponinC_Ca_Mg in component cytosolic_Ca_buffer (millimolar)"
    legend_rates[34] = "d/dt Mg_TroponinC_Ca_Mg in component cytosolic_Ca_buffer (millimolar)"
    legend_rates[35] = "d/dt Ca_Calmodulin in component cytosolic_Ca_buffer (millimolar)"
    legend_rates[36] = "d/dt Ca_Myosin in component cytosolic_Ca_buffer (millimolar)"
    legend_rates[37] = "d/dt Mg_Myosin in component cytosolic_Ca_buffer (millimolar)"
    legend_rates[38] = "d/dt Ca_SRB in component cytosolic_Ca_buffer (millimolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -8.556885e1
    constants[0] = 100
    constants[1] = 1000
    constants[2] = 5
    constants[3] = 9.5
    constants[4] = 5.4
    constants[5] = 140
    constants[6] = 1.8
    constants[7] = 150
    constants[8] = 135
    constants[9] = 1
    constants[10] = 15
    constants[11] = 8314.3
    constants[12] = 310
    constants[13] = 96485
    constants[14] = 1.381e-10
    constants[15] = 100
    constants[16] = 10.25
    states[1] = 8.80853
    states[2] = 8.80329
    states[3] = 8.80733
    states[4] = 8.597401e-5
    states[5] = 1.737475e-4
    states[6] = 1.031812e-4
    constants[17] = 16
    constants[18] = 0.11
    constants[19] = 0.89
    states[7] = 1.405627e-3
    states[8] = 9.867005e-1
    states[9] = 9.91562e-1
    constants[20] = 0.297e-3
    constants[21] = 0.11
    constants[22] = 0.89
    constants[23] = 4
    constants[24] = 11
    constants[25] = 1.5
    constants[26] = 1.90719
    constants[27] = 0.11
    constants[28] = 0.89
    states[10] = 8.641386e-3
    constants[29] = 0.11
    constants[30] = 0.89
    constants[31] = 0.01833
    states[11] = 5.412034e-3
    constants[32] = 0.001
    constants[33] = 0.06
    states[12] = 9.945511e-1
    states[13] = 4.051574e-3
    states[14] = 0.9946
    constants[34] = 0.02
    states[15] = 9.945511e-1
    states[16] = 4.051574e-3
    constants[35] = 0.109625
    constants[36] = 0.1
    constants[37] = 0.11
    constants[38] = 0.89
    constants[39] = 0.009
    constants[40] = 5.4e-4
    constants[41] = 1.5e-8
    constants[42] = 2.7e-7
    constants[43] = 0.9
    constants[44] = 0.1
    constants[45] = 0.341
    constants[46] = 0.341
    constants[47] = 0.75
    constants[48] = 0.75
    constants[49] = 0.75
    constants[50] = 0.75
    constants[51] = 1.8
    states[17] = 7.175662e-6
    states[18] = 1.000681
    states[19] = 1.452605e-2
    states[20] = 2.421991e-2
    constants[52] = 9
    constants[53] = 0.11
    constants[54] = 0.89
    constants[55] = 1.57
    constants[56] = 12.29
    constants[57] = 1.3
    constants[58] = 87.5
    constants[59] = 0.00359
    constants[60] = 0.000256
    constants[61] = 0.27
    constants[62] = 0.35
    constants[63] = 3
    constants[64] = 0.11
    constants[65] = 0.89
    constants[66] = 2.35
    constants[67] = 0.0005
    constants[68] = 1.6
    constants[69] = 0.0673
    constants[70] = 0.0002513
    constants[71] = 0.11
    constants[72] = 0.89
    constants[73] = 15
    constants[74] = 1
    constants[75] = 0.45
    states[21] = 8.884332e-1
    states[22] = 1.024274e-7
    states[23] = 8.156628e-7
    constants[76] = 25
    constants[77] = 10
    constants[78] = 0.06
    constants[79] = 0.5
    constants[80] = 0.005
    constants[81] = 2.5
    states[24] = 5.545201e-1
    constants[82] = 5.348e-6
    constants[83] = 5.3114e-3
    constants[84] = 2.6
    constants[85] = 0.000246
    constants[86] = 1.7
    constants[87] = 1.787
    states[25] = 7.720854e-1
    states[26] = 3.539892
    constants[88] = 1.65
    constants[89] = 7.561
    constants[90] = 0.0001
    constants[91] = 1e-3
    states[27] = 1.110363e-1
    states[28] = 9.566355e-3
    states[29] = 7.297378e-2
    states[30] = 7.347888e-3
    states[31] = 1.242988
    constants[92] = 0.0374
    constants[93] = 0.0046
    constants[94] = 0.0134
    constants[95] = 0.00165
    constants[96] = 0.14
    constants[97] = 100
    constants[98] = 100
    constants[99] = 1.3
    constants[100] = 30e-3
    constants[101] = 65
    states[32] = 8.773191e-3
    states[33] = 1.078283e-1
    states[34] = 1.524002e-2
    states[35] = 2.911916e-4
    states[36] = 1.298754e-3
    states[37] = 1.381982e-1
    states[38] = 2.143165e-3
    constants[102] = 0.07
    constants[103] = 0.14
    constants[104] = 0.14
    constants[105] = 0.024
    constants[106] = 0.14
    constants[107] = 0.14
    constants[108] = 0.0171
    constants[109] = 32.7
    constants[110] = 2.37
    constants[111] = 3e-3
    constants[112] = 34
    constants[113] = 13.8
    constants[114] = 15.7e-3
    constants[115] = 100
    constants[116] = 19.6e-3
    constants[117] = 0.032e-3
    constants[118] = 3.33e-3
    constants[119] = 238e-3
    constants[120] = 0.46e-3
    constants[121] = 0.057e-3
    constants[122] = 60e-3
    constants[123] = ((constants[11]*constants[12])/constants[13])*log(constants[4]/constants[8])
    constants[124] = ((constants[11]*constants[12])/constants[13])*log(constants[10]/constants[7])
    constants[125] = (exp(constants[5]/67.3000)-1.00000)/7.00000
    constants[126] = 0.0300000*(power(constants[4]/5.40000, 1.0/2))
    constants[127] = 0.900000*(power(constants[4]/5.40000, 1.0/2))
    constants[128] = power(constants[51], (constants[12]-310.000)/10.0000)
    constants[129] = power(constants[55], (constants[12]-310.000)/10.0000)
    constants[130] = power(constants[66], (constants[12]-310.000)/10.0000)
    constants[131] = power(constants[84], (constants[12]-310.000)/10.0000)
    constants[132] = (3.14159*(power(constants[16]/1000.00, 2.00000))*constants[15])/(power(1000.00, 3.00000))
    constants[133] = 0.0350000*constants[132]
    constants[134] = 0.0200000*constants[132]
    constants[135] = 0.0539000*0.0100000*constants[132]
    constants[136] = 0.650000*constants[132]
    return (states, constants)

def computeRates(voi, states, constants,activation): #activation
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[19] = 1.70000*states[6]*(1.00000-states[19])-0.0119000*states[19]
    rates[20] = 1.70000*states[5]*(1.00000-states[20])-0.0119000*states[20]
    algebraic[0] = custom_piecewise([less(states[0] , -40.0000), 0.135000*exp((80.0000+states[0])/-6.80000) , True, 0.00000])
    algebraic[14] = custom_piecewise([less(states[0] , -40.0000), 3.56000*exp(0.0790000*states[0])+310000.*exp(0.350000*states[0]) , True, 1.00000/(0.130000*(1.00000+exp((states[0]+10.6600)/-11.1000)))])
    rates[8] = algebraic[0]*(1.00000-states[8])-algebraic[14]*states[8]
    algebraic[1] = custom_piecewise([less(states[0] , -40.0000), (((-127140.*exp(0.244400*states[0])-3.47400e-05*exp(-0.0439100*states[0]))*(states[0]+37.7800))/1.00000)/(1.00000+exp(0.311000*(states[0]+79.2300))) , True, 0.00000])
    algebraic[15] = custom_piecewise([less(states[0] , -40.0000), (0.121200*exp(-0.0105200*states[0]))/(1.00000+exp(-0.137800*(states[0]+40.1400))) , True, (0.300000*exp(-2.53500e-07*states[0]))/(1.00000+exp(-0.100000*(states[0]+32.0000)))])
    rates[9] = algebraic[1]*(1.00000-states[9])-algebraic[15]*states[9]
    algebraic[2] = ((0.320000*(states[0]+47.1300))/1.00000)/(1.00000-exp(-0.100000*(states[0]+47.1300)))
    algebraic[16] = 0.0800000*exp(-states[0]/11.0000)
    rates[7] = algebraic[2]*(1.00000-states[7])-algebraic[16]*states[7]
    algebraic[3] = 1.00000/(1.00000+exp(-(50.0000+states[0])/7.50000))
    algebraic[17] = 1.00000/((0.00138000*(states[0]+7.00000))/(1.00000-exp(-0.123000*(states[0]+7.00000)))+(0.000610000*(states[0]+10.0000))/(exp(0.145000*(states[0]+10.0000))-1.00000))
    rates[10] = (algebraic[3]-states[10])/algebraic[17]
    algebraic[5] = 1.00000/(1.00000+exp(-(states[0]-1.50000)/16.7000))
    algebraic[19] = 1.00000/((7.19000e-05*(states[0]+30.0000))/(1.00000-exp(-0.148000*(states[0]+30.0000)))+(0.000131000*(states[0]+30.0000))/(-1.00000+exp(0.0687000*(states[0]+30.0000))))
    rates[11] = (algebraic[5]-states[11])/algebraic[19]
    algebraic[6] = 1.00000/(1.00000+exp(-(states[0]+3.00000)/15.0000))
    algebraic[20] = 9.00000/(1.00000+exp((states[0]+3.00000)/15.0000))+0.500000
    rates[13] = (algebraic[6]-states[13])/algebraic[20]
    algebraic[7] = 1.00000/(1.00000+exp((states[0]+33.5000)/10.0000))
    algebraic[21] = 3000.00/(1.00000+exp((states[0]+60.0000)/10.0000))+30.0000
    rates[12] = (algebraic[7]-states[12])/algebraic[21]
    algebraic[8] = 1.00000/(1.00000+exp((states[0]+33.5000)/10.0000))
    algebraic[22] = 2800.00/(1.00000+exp((states[0]+60.0000)/10.0000))+220.000
    rates[14] = (algebraic[8]-states[14])/algebraic[22]
    algebraic[9] = 1.00000/(1.00000+exp(-(states[0]+3.00000)/15.0000))
    algebraic[23] = 3.50000*exp(-(power(states[0]/30.0000, 2.00000)))+1.50000
    rates[16] = (algebraic[9]-states[16])/algebraic[23]
    algebraic[10] = 1.00000/(1.00000+exp((states[0]+33.5000)/10.0000))
    algebraic[24] = 20.0000/(1.00000+exp((states[0]+33.5000)/10.0000))+20.0000
    rates[15] = (algebraic[10]-states[15])/algebraic[24]
    algebraic[11] = 1.00000/(1.00000+exp(-(states[0]+14.5000)/6.00000))
    algebraic[25] = (1.00000*algebraic[11]*(1.00000-exp(-(states[0]+14.5000)/6.00000)))/(0.0350000*(states[0]+14.5000))
    rates[17] = (algebraic[11]-states[17])/algebraic[25]
    algebraic[12] = 1.00000/(1.00000+exp((states[0]+35.0600)/3.60000))+0.600000/(1.00000+exp((50.0000-states[0])/20.0000))
    algebraic[26] = 1.00000/(0.0197000*exp(-(power(0.0337000*(states[0]+14.5000), 2.00000)))+0.0200000)
    rates[18] = (algebraic[12]-states[18])/algebraic[26]
    algebraic[13] = constants[73]-(constants[73]-constants[74])/(1.00000+power(constants[75]/states[24], constants[81]))
    algebraic[27] = constants[77]/algebraic[13]
    algebraic[29] = constants[79]*algebraic[13]
    rates[23] = (algebraic[27]*(power(states[5], 2.00000))*states[21]-constants[78]*states[23])-(algebraic[29]*states[5]*states[23]-constants[80]*states[22])
    algebraic[31] = ((1.00000-states[21])-states[23])-states[22]
    rates[21] = (constants[80]*algebraic[31]-algebraic[29]*states[5]*states[21])-(algebraic[27]*(power(states[5], 2.00000))*states[21]-constants[78]*states[23])
    rates[22] = (algebraic[29]*states[5]*states[23]-constants[80]*states[22])-(constants[78]*states[22]-algebraic[27]*(power(states[5], 2.00000))*algebraic[31])
    algebraic[80] = constants[90]*states[2]*(constants[89]-states[26])-constants[91]*states[26]
    rates[26] = algebraic[80]
    algebraic[18] = ((constants[11]*constants[12])/constants[13])*log(constants[5]/states[2])
    algebraic[33] = (power(states[7], 3.00000))*states[8]*states[9]
    algebraic[34] = constants[18]*constants[17]*algebraic[33]*(states[0]-algebraic[18])
    algebraic[37] = constants[21]*constants[20]*(states[0]-algebraic[18])
    algebraic[40] = 1.00000/(1.00000+0.124500*exp((-0.100000*states[0]*constants[13])/(constants[11]*constants[12]))+0.0365000*constants[125]*exp((-states[0]*constants[13])/(constants[11]*constants[12])))
    algebraic[41] = (((constants[27]*constants[26]*algebraic[40])/(1.00000+power(constants[24]/states[2], constants[23])))*constants[4])/(constants[4]+constants[25])
    algebraic[67] = 1.00000-states[20]
    algebraic[63] = (0.450000*states[17]*states[18]*constants[128]*states[0]*(power(constants[13], 2.00000)))/(constants[11]*constants[12])
    algebraic[69] = (algebraic[63]*algebraic[67]*constants[43]*constants[41]*(constants[47]*states[2]*exp((states[0]*constants[13])/(constants[11]*constants[12]))-constants[48]*constants[5]))/(exp((states[0]*constants[13])/(constants[11]*constants[12]))-1.00000)
    algebraic[75] = 1.00000/(1.00000+power(constants[60]/states[5], 3.00000))
    algebraic[72] = (exp((constants[62]*states[0]*constants[13])/(constants[11]*constants[12]))*(power(states[2], constants[63]))*constants[6]-exp(((constants[62]-1.00000)*states[0]*constants[13])/(constants[11]*constants[12]))*(power(constants[5], constants[63]))*states[5])/(1.00000+constants[61]*exp(((constants[62]-1.00000)*states[0]*constants[13])/(constants[11]*constants[12])))
    algebraic[76] = (constants[53]*constants[52]*algebraic[75]*constants[129]*algebraic[72])/(constants[59]*(power(constants[5], constants[63]))*(1.00000+power(states[2]/constants[56], constants[63]))+(power(constants[58], constants[63]))*states[5]*(1.00000+states[5]/constants[59])+constants[57]*(power(states[2], constants[63]))+(power(states[2], constants[63]))*constants[6]+(power(constants[5], constants[63]))*states[5])
    algebraic[78] = (states[2]-states[3])*1.83130e-14
    rates[2] = ((-constants[14]*(algebraic[34]+3.00000*algebraic[76]+algebraic[37]+3.00000*algebraic[41]+algebraic[69]))/(constants[135]*constants[13])-algebraic[78]/constants[135])-algebraic[80]
    algebraic[81] = (states[3]-states[1])*1.63860e-12
    rates[1] = algebraic[81]/constants[136]
    algebraic[83] = constants[90]*states[3]*(constants[88]-states[25])-constants[91]*states[25]
    rates[25] = algebraic[83]
    algebraic[28] = ((constants[11]*constants[12])/constants[13])*log(constants[5]/states[3])
    algebraic[35] = constants[19]*constants[17]*algebraic[33]*(states[0]-algebraic[28])
    algebraic[38] = constants[22]*constants[20]*(states[0]-algebraic[28])
    algebraic[42] = (((constants[28]*constants[26]*algebraic[40])/(1.00000+power(constants[24]/states[3], constants[23])))*constants[4])/(constants[4]+constants[25])
    algebraic[64] = 1.00000-states[19]
    algebraic[66] = (algebraic[63]*algebraic[64]*constants[44]*constants[41]*(constants[47]*states[3]*exp((states[0]*constants[13])/(constants[11]*constants[12]))-constants[48]*constants[5]))/(exp((states[0]*constants[13])/(constants[11]*constants[12]))-1.00000)
    algebraic[74] = 1.00000/(1.00000+power(constants[60]/states[6], 3.00000))
    algebraic[73] = (exp((constants[62]*states[0]*constants[13])/(constants[11]*constants[12]))*(power(states[3], constants[63]))*constants[6]-exp(((constants[62]-1.00000)*states[0]*constants[13])/(constants[11]*constants[12]))*(power(constants[5], constants[63]))*states[6])/(1.00000+constants[61]*exp(((constants[62]-1.00000)*states[0]*constants[13])/(constants[11]*constants[12])))
    algebraic[77] = (constants[54]*constants[52]*algebraic[74]*constants[129]*algebraic[73])/(constants[59]*(power(constants[5], constants[63]))*(1.00000+power(states[3]/constants[56], constants[63]))+(power(constants[58], constants[63]))*states[6]*(1.00000+states[6]/constants[59])+constants[57]*(power(states[3], constants[63]))+(power(states[3], constants[63]))*constants[6]+(power(constants[5], constants[63]))*states[6])
    rates[3] = ((-constants[14]*(algebraic[35]+3.00000*algebraic[77]+algebraic[38]+3.00000*algebraic[42]+algebraic[66]))/(constants[134]*constants[13])+(algebraic[78]-algebraic[81])/constants[134])-algebraic[83]
    algebraic[36] = algebraic[34]+algebraic[35]
    algebraic[39] = algebraic[37]+algebraic[38]
    algebraic[43] = algebraic[41]+algebraic[42]
    algebraic[44] = 1.00000/(1.00000+exp((33.0000+states[0])/22.4000))
    algebraic[45] = constants[126]*states[10]*algebraic[44]*(states[0]-constants[123])
    algebraic[46] = -log10(states[5]/1.00000)+3.00000
    algebraic[48] = 0.0700000*(0.0570000+0.190000/(1.00000+exp((-7.20000+algebraic[46])/0.600000)))
    algebraic[50] = ((constants[11]*constants[12])/constants[13])*log((constants[4]+constants[31]*constants[5])/(constants[8]+constants[31]*states[1]))
    algebraic[51] = constants[29]*algebraic[48]*(power(states[11], 2.00000))*(states[0]-algebraic[50])
    algebraic[47] = -log10(states[6]/1.00000)+3.00000
    algebraic[49] = 0.0700000*(0.0570000+0.190000/(1.00000+exp((-7.20000+algebraic[47])/0.600000)))
    algebraic[52] = constants[30]*algebraic[49]*(power(states[11], 2.00000))*(states[0]-algebraic[50])
    algebraic[53] = algebraic[51]+algebraic[52]
    algebraic[54] = (constants[32]*(states[0]-constants[123]))/(1.00000+exp(7.48800-states[0]/5.98000))
    algebraic[55] = constants[33]*states[13]*(states[12]+0.500000*states[14])*(states[0]-constants[123])
    algebraic[56] = constants[34]*states[16]*states[15]*(states[0]-constants[123])
    algebraic[57] = 1.02000/(1.00000+exp(0.238500*(states[0]-(constants[123]+59.2150))))
    algebraic[58] = (0.491240*exp(0.0803200*((states[0]-constants[123])+5.47600))+1.00000*exp(0.0617500*(states[0]-(constants[123]+594.310))))/(1.00000+exp(-0.514300*((states[0]-constants[123])+4.75300)))
    algebraic[59] = algebraic[57]/(algebraic[57]+algebraic[58])
    algebraic[60] = constants[127]*algebraic[59]*(states[0]-constants[123])
    algebraic[79] = algebraic[76]+algebraic[77]
    algebraic[61] = constants[35]*(states[0]-constants[124])*(constants[37]/(1.00000+constants[36]/states[5])+constants[38]/(1.00000+constants[36]/states[6]))
    algebraic[62] = constants[39]*(states[0]-constants[124])
    algebraic[65] = (algebraic[63]*algebraic[64]*constants[44]*constants[40]*4.00000*(constants[45]*states[6]*exp((2.00000*states[0]*constants[13])/(constants[11]*constants[12]))-constants[46]*constants[6]))/(exp((2.00000*states[0]*constants[13])/(constants[11]*constants[12]))-1.00000)
    algebraic[68] = (algebraic[63]*algebraic[67]*constants[43]*constants[40]*4.00000*(constants[45]*states[5]*exp((2.00000*states[0]*constants[13])/(constants[11]*constants[12]))-constants[46]*constants[6]))/(exp((2.00000*states[0]*constants[13])/(constants[11]*constants[12]))-1.00000)
    algebraic[70] = (algebraic[63]*(algebraic[64]*constants[44]+algebraic[67]*constants[43])*constants[42]*(constants[49]*constants[8]*exp((states[0]*constants[13])/(constants[11]*constants[12]))-constants[50]*constants[4]))/(exp((states[0]*constants[13])/(constants[11]*constants[12]))-1.00000)
    algebraic[71] = algebraic[65]+algebraic[68]+algebraic[66]+algebraic[69]+algebraic[70]
    algebraic[30] = ((constants[11]*constants[12])/(2.00000*constants[13]))*log(constants[6]/states[5])
    algebraic[86] = constants[70]*constants[71]*(states[0]-algebraic[30])
    algebraic[32] = ((constants[11]*constants[12])/(2.00000*constants[13]))*log(constants[6]/states[6])
    algebraic[87] = constants[70]*constants[72]*(states[0]-algebraic[32])
    algebraic[88] = algebraic[87]+algebraic[86]
    algebraic[84] = (constants[130]*constants[69]*constants[65])/(1.00000+power(constants[67]/states[6], constants[68]))
    algebraic[82] = (constants[130]*constants[69]*constants[64])/(1.00000+power(constants[67]/states[5], constants[68]))
    algebraic[85] = algebraic[82]+algebraic[84]

    #print('voi=',voi)
    algebraic[4] = custom_piecewise([greater_equal(voi-floor(voi/constants[1])*constants[1] , \
                   constants[0]) & less_equal(voi-floor(voi/constants[1])*constants[1] , \
                   constants[0]+constants[2]), -constants[3] , True, 0.00000])
    print('algebraic[4]',algebraic[4])
    #algebraic[4] = activation * -9.5
    rates[0] = -(algebraic[36]+algebraic[39]+algebraic[43]+algebraic[45]\
                +algebraic[53]+algebraic[55]+algebraic[56]+algebraic[60]+algebraic[79]\
                +algebraic[61]+algebraic[62]+algebraic[71]+algebraic[88]+algebraic[85]+algebraic[54]+algebraic[4])
    algebraic[92] = constants[98]*states[24]*((constants[96]*constants[136])/constants[133]-states[31])-constants[101]*states[31]
    rates[31] = algebraic[92]
    algebraic[89] = constants[76]*states[23]*(states[24]-states[5])
    algebraic[90] = constants[82]*(states[24]-states[5])
    algebraic[91] = (constants[131]*constants[83]*(power(states[4]/constants[85], constants[87])-power(states[24]/constants[86], constants[87])))/(1.00000+power(states[4]/constants[85], constants[87])+power(states[24]/constants[86], constants[87]))
    rates[24] = (algebraic[91]-((algebraic[90]*constants[136])/constants[133]+algebraic[89]))-algebraic[92]
    algebraic[95] = constants[97]*states[5]*((constants[93]*0.100000*constants[136])/constants[135]-states[28])-constants[99]*states[28]
    rates[28] = algebraic[95]
    algebraic[97] = constants[97]*states[6]*((constants[92]*constants[136])/constants[134]-states[27])-constants[99]*states[27]
    rates[27] = algebraic[97]
    algebraic[98] = constants[97]*states[5]*((constants[95]*0.100000*constants[136])/constants[135]-states[30])-constants[100]*states[30]
    rates[30] = algebraic[98]
    algebraic[96] = constants[109]*states[4]*(constants[102]-states[32])-constants[116]*states[32]
    rates[32] = algebraic[96]
    algebraic[100] = constants[97]*states[6]*((constants[94]*constants[136])/constants[134]-states[29])-constants[100]*states[29]
    rates[29] = algebraic[100]
    algebraic[99] = constants[110]*states[4]*(constants[103]-(states[33]+states[34]))-constants[117]*states[33]
    rates[33] = algebraic[99]
    algebraic[93] = (states[5]-states[6])*8.24130e-13
    algebraic[101] = algebraic[95]+algebraic[98]
    algebraic[103] = (algebraic[68]-2.00000*algebraic[76])+algebraic[86]+algebraic[82]
    rates[5] = (((-algebraic[103]*constants[14])/(constants[135]*2.00000*constants[13])-algebraic[93]/constants[135])+(algebraic[89]*constants[133])/constants[135]+(algebraic[90]*constants[136])/constants[135])-1.00000*algebraic[101]
    algebraic[102] = constants[111]*constants[9]*(constants[104]-(states[33]+states[34]))-constants[118]*states[34]
    rates[34] = algebraic[102]
    algebraic[94] = (states[6]-states[4])*3.72430e-12
    algebraic[104] = algebraic[97]+algebraic[100]
    algebraic[106] = (algebraic[65]-2.00000*algebraic[77])+algebraic[87]+algebraic[84]
    rates[6] = ((-algebraic[106]*constants[14])/(constants[134]*2.00000*constants[13])+(algebraic[93]-algebraic[94])/constants[134])-1.00000*algebraic[104]
    algebraic[105] = constants[112]*states[4]*(constants[105]-states[35])-constants[119]*states[35]
    rates[35] = algebraic[105]
    algebraic[107] = constants[113]*states[4]*(constants[106]-(states[36]+states[37]))-constants[120]*states[36]
    rates[36] = algebraic[107]
    algebraic[108] = constants[114]*constants[9]*(constants[107]-(states[36]+states[37]))-constants[121]*states[37]
    rates[37] = algebraic[108]
    algebraic[109] = constants[115]*states[4]*(constants[108]-states[38])-constants[122]*states[38]
    rates[38] = algebraic[109]
    algebraic[110] = algebraic[96]+algebraic[99]+algebraic[102]+algebraic[105]+algebraic[107]+algebraic[108]+algebraic[109]
    rates[4] = ((-algebraic[91]*constants[133])/constants[136]+algebraic[94]/constants[136])-1.00000*algebraic[110]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = custom_piecewise([less(states[0] , -40.0000), 0.135000*exp((80.0000+states[0])/-6.80000) , True, 0.00000])
    algebraic[14] = custom_piecewise([less(states[0] , -40.0000), 3.56000*exp(0.0790000*states[0])+310000.*exp(0.350000*states[0]) , True, 1.00000/(0.130000*(1.00000+exp((states[0]+10.6600)/-11.1000)))])
    algebraic[1] = custom_piecewise([less(states[0] , -40.0000), (((-127140.*exp(0.244400*states[0])-3.47400e-05*exp(-0.0439100*states[0]))*(states[0]+37.7800))/1.00000)/(1.00000+exp(0.311000*(states[0]+79.2300))) , True, 0.00000])
    algebraic[15] = custom_piecewise([less(states[0] , -40.0000), (0.121200*exp(-0.0105200*states[0]))/(1.00000+exp(-0.137800*(states[0]+40.1400))) , True, (0.300000*exp(-2.53500e-07*states[0]))/(1.00000+exp(-0.100000*(states[0]+32.0000)))])
    algebraic[2] = ((0.320000*(states[0]+47.1300))/1.00000)/(1.00000-exp(-0.100000*(states[0]+47.1300)))
    algebraic[16] = 0.0800000*exp(-states[0]/11.0000)
    algebraic[3] = 1.00000/(1.00000+exp(-(50.0000+states[0])/7.50000))
    algebraic[17] = 1.00000/((0.00138000*(states[0]+7.00000))/(1.00000-exp(-0.123000*(states[0]+7.00000)))+(0.000610000*(states[0]+10.0000))/(exp(0.145000*(states[0]+10.0000))-1.00000))
    algebraic[5] = 1.00000/(1.00000+exp(-(states[0]-1.50000)/16.7000))
    algebraic[19] = 1.00000/((7.19000e-05*(states[0]+30.0000))/(1.00000-exp(-0.148000*(states[0]+30.0000)))+(0.000131000*(states[0]+30.0000))/(-1.00000+exp(0.0687000*(states[0]+30.0000))))
    algebraic[6] = 1.00000/(1.00000+exp(-(states[0]+3.00000)/15.0000))
    algebraic[20] = 9.00000/(1.00000+exp((states[0]+3.00000)/15.0000))+0.500000
    algebraic[7] = 1.00000/(1.00000+exp((states[0]+33.5000)/10.0000))
    algebraic[21] = 3000.00/(1.00000+exp((states[0]+60.0000)/10.0000))+30.0000
    algebraic[8] = 1.00000/(1.00000+exp((states[0]+33.5000)/10.0000))
    algebraic[22] = 2800.00/(1.00000+exp((states[0]+60.0000)/10.0000))+220.000
    algebraic[9] = 1.00000/(1.00000+exp(-(states[0]+3.00000)/15.0000))
    algebraic[23] = 3.50000*exp(-(power(states[0]/30.0000, 2.00000)))+1.50000
    algebraic[10] = 1.00000/(1.00000+exp((states[0]+33.5000)/10.0000))
    algebraic[24] = 20.0000/(1.00000+exp((states[0]+33.5000)/10.0000))+20.0000
    algebraic[11] = 1.00000/(1.00000+exp(-(states[0]+14.5000)/6.00000))
    algebraic[25] = (1.00000*algebraic[11]*(1.00000-exp(-(states[0]+14.5000)/6.00000)))/(0.0350000*(states[0]+14.5000))
    algebraic[12] = 1.00000/(1.00000+exp((states[0]+35.0600)/3.60000))+0.600000/(1.00000+exp((50.0000-states[0])/20.0000))
    algebraic[26] = 1.00000/(0.0197000*exp(-(power(0.0337000*(states[0]+14.5000), 2.00000)))+0.0200000)
    algebraic[13] = constants[73]-(constants[73]-constants[74])/(1.00000+power(constants[75]/states[24], constants[81]))
    algebraic[27] = constants[77]/algebraic[13]
    algebraic[29] = constants[79]*algebraic[13]
    algebraic[31] = ((1.00000-states[21])-states[23])-states[22]
    algebraic[80] = constants[90]*states[2]*(constants[89]-states[26])-constants[91]*states[26]
    algebraic[18] = ((constants[11]*constants[12])/constants[13])*log(constants[5]/states[2])
    algebraic[33] = (power(states[7], 3.00000))*states[8]*states[9]
    algebraic[34] = constants[18]*constants[17]*algebraic[33]*(states[0]-algebraic[18])
    algebraic[37] = constants[21]*constants[20]*(states[0]-algebraic[18])
    algebraic[40] = 1.00000/(1.00000+0.124500*exp((-0.100000*states[0]*constants[13])/(constants[11]*constants[12]))+0.0365000*constants[125]*exp((-states[0]*constants[13])/(constants[11]*constants[12])))
    algebraic[41] = (((constants[27]*constants[26]*algebraic[40])/(1.00000+power(constants[24]/states[2], constants[23])))*constants[4])/(constants[4]+constants[25])
    algebraic[67] = 1.00000-states[20]
    algebraic[63] = (0.450000*states[17]*states[18]*constants[128]*states[0]*(power(constants[13], 2.00000)))/(constants[11]*constants[12])
    algebraic[69] = (algebraic[63]*algebraic[67]*constants[43]*constants[41]*(constants[47]*states[2]*exp((states[0]*constants[13])/(constants[11]*constants[12]))-constants[48]*constants[5]))/(exp((states[0]*constants[13])/(constants[11]*constants[12]))-1.00000)
    algebraic[75] = 1.00000/(1.00000+power(constants[60]/states[5], 3.00000))
    algebraic[72] = (exp((constants[62]*states[0]*constants[13])/(constants[11]*constants[12]))*(power(states[2], constants[63]))*constants[6]-exp(((constants[62]-1.00000)*states[0]*constants[13])/(constants[11]*constants[12]))*(power(constants[5], constants[63]))*states[5])/(1.00000+constants[61]*exp(((constants[62]-1.00000)*states[0]*constants[13])/(constants[11]*constants[12])))
    algebraic[76] = (constants[53]*constants[52]*algebraic[75]*constants[129]*algebraic[72])/(constants[59]*(power(constants[5], constants[63]))*(1.00000+power(states[2]/constants[56], constants[63]))+(power(constants[58], constants[63]))*states[5]*(1.00000+states[5]/constants[59])+constants[57]*(power(states[2], constants[63]))+(power(states[2], constants[63]))*constants[6]+(power(constants[5], constants[63]))*states[5])
    algebraic[78] = (states[2]-states[3])*1.83130e-14
    algebraic[81] = (states[3]-states[1])*1.63860e-12
    algebraic[83] = constants[90]*states[3]*(constants[88]-states[25])-constants[91]*states[25]
    algebraic[28] = ((constants[11]*constants[12])/constants[13])*log(constants[5]/states[3])
    algebraic[35] = constants[19]*constants[17]*algebraic[33]*(states[0]-algebraic[28])
    algebraic[38] = constants[22]*constants[20]*(states[0]-algebraic[28])
    algebraic[42] = (((constants[28]*constants[26]*algebraic[40])/(1.00000+power(constants[24]/states[3], constants[23])))*constants[4])/(constants[4]+constants[25])
    algebraic[64] = 1.00000-states[19]
    algebraic[66] = (algebraic[63]*algebraic[64]*constants[44]*constants[41]*(constants[47]*states[3]*exp((states[0]*constants[13])/(constants[11]*constants[12]))-constants[48]*constants[5]))/(exp((states[0]*constants[13])/(constants[11]*constants[12]))-1.00000)
    algebraic[74] = 1.00000/(1.00000+power(constants[60]/states[6], 3.00000))
    algebraic[73] = (exp((constants[62]*states[0]*constants[13])/(constants[11]*constants[12]))*(power(states[3], constants[63]))*constants[6]-exp(((constants[62]-1.00000)*states[0]*constants[13])/(constants[11]*constants[12]))*(power(constants[5], constants[63]))*states[6])/(1.00000+constants[61]*exp(((constants[62]-1.00000)*states[0]*constants[13])/(constants[11]*constants[12])))
    algebraic[77] = (constants[54]*constants[52]*algebraic[74]*constants[129]*algebraic[73])/(constants[59]*(power(constants[5], constants[63]))*(1.00000+power(states[3]/constants[56], constants[63]))+(power(constants[58], constants[63]))*states[6]*(1.00000+states[6]/constants[59])+constants[57]*(power(states[3], constants[63]))+(power(states[3], constants[63]))*constants[6]+(power(constants[5], constants[63]))*states[6])
    algebraic[36] = algebraic[34]+algebraic[35]
    algebraic[39] = algebraic[37]+algebraic[38]
    algebraic[43] = algebraic[41]+algebraic[42]
    algebraic[44] = 1.00000/(1.00000+exp((33.0000+states[0])/22.4000))
    algebraic[45] = constants[126]*states[10]*algebraic[44]*(states[0]-constants[123])
    algebraic[46] = -log10(states[5]/1.00000)+3.00000
    algebraic[48] = 0.0700000*(0.0570000+0.190000/(1.00000+exp((-7.20000+algebraic[46])/0.600000)))
    algebraic[50] = ((constants[11]*constants[12])/constants[13])*log((constants[4]+constants[31]*constants[5])/(constants[8]+constants[31]*states[1]))
    algebraic[51] = constants[29]*algebraic[48]*(power(states[11], 2.00000))*(states[0]-algebraic[50])
    algebraic[47] = -log10(states[6]/1.00000)+3.00000
    algebraic[49] = 0.0700000*(0.0570000+0.190000/(1.00000+exp((-7.20000+algebraic[47])/0.600000)))
    algebraic[52] = constants[30]*algebraic[49]*(power(states[11], 2.00000))*(states[0]-algebraic[50])
    algebraic[53] = algebraic[51]+algebraic[52]
    algebraic[54] = (constants[32]*(states[0]-constants[123]))/(1.00000+exp(7.48800-states[0]/5.98000))
    algebraic[55] = constants[33]*states[13]*(states[12]+0.500000*states[14])*(states[0]-constants[123])
    algebraic[56] = constants[34]*states[16]*states[15]*(states[0]-constants[123])
    algebraic[57] = 1.02000/(1.00000+exp(0.238500*(states[0]-(constants[123]+59.2150))))
    algebraic[58] = (0.491240*exp(0.0803200*((states[0]-constants[123])+5.47600))+1.00000*exp(0.0617500*(states[0]-(constants[123]+594.310))))/(1.00000+exp(-0.514300*((states[0]-constants[123])+4.75300)))
    algebraic[59] = algebraic[57]/(algebraic[57]+algebraic[58])
    algebraic[60] = constants[127]*algebraic[59]*(states[0]-constants[123])
    algebraic[79] = algebraic[76]+algebraic[77]
    algebraic[61] = constants[35]*(states[0]-constants[124])*(constants[37]/(1.00000+constants[36]/states[5])+constants[38]/(1.00000+constants[36]/states[6]))
    algebraic[62] = constants[39]*(states[0]-constants[124])
    algebraic[65] = (algebraic[63]*algebraic[64]*constants[44]*constants[40]*4.00000*(constants[45]*states[6]*exp((2.00000*states[0]*constants[13])/(constants[11]*constants[12]))-constants[46]*constants[6]))/(exp((2.00000*states[0]*constants[13])/(constants[11]*constants[12]))-1.00000)
    algebraic[68] = (algebraic[63]*algebraic[67]*constants[43]*constants[40]*4.00000*(constants[45]*states[5]*exp((2.00000*states[0]*constants[13])/(constants[11]*constants[12]))-constants[46]*constants[6]))/(exp((2.00000*states[0]*constants[13])/(constants[11]*constants[12]))-1.00000)
    algebraic[70] = (algebraic[63]*(algebraic[64]*constants[44]+algebraic[67]*constants[43])*constants[42]*(constants[49]*constants[8]*exp((states[0]*constants[13])/(constants[11]*constants[12]))-constants[50]*constants[4]))/(exp((states[0]*constants[13])/(constants[11]*constants[12]))-1.00000)
    algebraic[71] = algebraic[65]+algebraic[68]+algebraic[66]+algebraic[69]+algebraic[70]
    algebraic[30] = ((constants[11]*constants[12])/(2.00000*constants[13]))*log(constants[6]/states[5])
    algebraic[86] = constants[70]*constants[71]*(states[0]-algebraic[30])
    algebraic[32] = ((constants[11]*constants[12])/(2.00000*constants[13]))*log(constants[6]/states[6])
    algebraic[87] = constants[70]*constants[72]*(states[0]-algebraic[32])
    algebraic[88] = algebraic[87]+algebraic[86]
    algebraic[84] = (constants[130]*constants[69]*constants[65])/(1.00000+power(constants[67]/states[6], constants[68]))
    algebraic[82] = (constants[130]*constants[69]*constants[64])/(1.00000+power(constants[67]/states[5], constants[68]))
    algebraic[85] = algebraic[82]+algebraic[84]
    algebraic[4] = custom_piecewise([greater_equal(voi-floor(voi/constants[1])*constants[1] , constants[0]) & less_equal(voi-floor(voi/constants[1])*constants[1] , constants[0]+constants[2]), -constants[3] , True, 0.00000])
    algebraic[92] = constants[98]*states[24]*((constants[96]*constants[136])/constants[133]-states[31])-constants[101]*states[31]
    algebraic[89] = constants[76]*states[23]*(states[24]-states[5])
    algebraic[90] = constants[82]*(states[24]-states[5])
    algebraic[91] = (constants[131]*constants[83]*(power(states[4]/constants[85], constants[87])-power(states[24]/constants[86], constants[87])))/(1.00000+power(states[4]/constants[85], constants[87])+power(states[24]/constants[86], constants[87]))
    algebraic[95] = constants[97]*states[5]*((constants[93]*0.100000*constants[136])/constants[135]-states[28])-constants[99]*states[28]
    algebraic[97] = constants[97]*states[6]*((constants[92]*constants[136])/constants[134]-states[27])-constants[99]*states[27]
    algebraic[98] = constants[97]*states[5]*((constants[95]*0.100000*constants[136])/constants[135]-states[30])-constants[100]*states[30]
    algebraic[96] = constants[109]*states[4]*(constants[102]-states[32])-constants[116]*states[32]
    algebraic[100] = constants[97]*states[6]*((constants[94]*constants[136])/constants[134]-states[29])-constants[100]*states[29]
    algebraic[99] = constants[110]*states[4]*(constants[103]-(states[33]+states[34]))-constants[117]*states[33]
    algebraic[93] = (states[5]-states[6])*8.24130e-13
    algebraic[101] = algebraic[95]+algebraic[98]
    algebraic[103] = (algebraic[68]-2.00000*algebraic[76])+algebraic[86]+algebraic[82]
    algebraic[102] = constants[111]*constants[9]*(constants[104]-(states[33]+states[34]))-constants[118]*states[34]
    algebraic[94] = (states[6]-states[4])*3.72430e-12
    algebraic[104] = algebraic[97]+algebraic[100]
    algebraic[106] = (algebraic[65]-2.00000*algebraic[77])+algebraic[87]+algebraic[84]
    algebraic[105] = constants[112]*states[4]*(constants[105]-states[35])-constants[119]*states[35]
    algebraic[107] = constants[113]*states[4]*(constants[106]-(states[36]+states[37]))-constants[120]*states[36]
    algebraic[108] = constants[114]*constants[9]*(constants[107]-(states[36]+states[37]))-constants[121]*states[37]
    algebraic[109] = constants[115]*states[4]*(constants[108]-states[38])-constants[122]*states[38]
    algebraic[110] = algebraic[96]+algebraic[99]+algebraic[102]+algebraic[105]+algebraic[107]+algebraic[108]+algebraic[109]
    return algebraic

def custom_piecewise(cases):
    """Compute result of a piecewise function"""
    return select(cases[0::2],cases[1::2])

def solve_model():
    """Solve model with ODE solver"""
    import numpy as np
    from scipy.integrate import ode
    from scipy.integrate import solve_ivp
    from functools import partial
    # Initialise constants and state variables
    (init_states, constants) = initConsts()

    #print("states=",init_states)
    #print('states size=',np.shape(init_states))
    #print('cosntants size=',np.shape(constants))

    # Set timespan to solve over
    voi = linspace(0, 10000, 10000)

    # Construct ODE object to solve
    r = ode(computeRates)
    r.set_integrator('vode', method='bdf', atol=1e-06, rtol=1e-06, max_step=1)
    r.set_initial_value(init_states, voi[0])
    r.set_f_params(constants)

    # Solve model
    states = array([[0.0] * len(voi)] * sizeStates)
    states[:,0] = init_states
    for (i,t) in enumerate(voi[1:]):
        if r.successful():
            r.integrate(t)
            states[:,i+1] = r.y
        else:
            break

    # Compute algebraic variables
    algebraic = computeAlgebraic(constants, states, voi)
    print("dimension states=",np.shape(states))
    print("dimension algebraic=",np.shape(algebraic))

    return (voi, states, algebraic)

def solve_system():
    from scipy.integrate import ode
    from scipy.integrate import solve_ivp
    import numpy as np
    from functools import partial

    (init_states, constants) = initConsts()

    # Set timespan to solve over
    dt=1
    voi = np.arange(0,1000,dt)

    states = array([[0.0] * len(voi)] * sizeStates)
    states[:,0] = init_states

    for (i,t) in enumerate(voi[1:]):
        if ((i>100)&(i<110)):
            activation = 1.0
        else:
            activation = 0.0
        sol = solve_ivp(partial(computeRates,
                            constants=constants,
                            activation=activation),\
                            voi[i]+[0, dt], states[:,i],
                    method='BDF')


        if (i<(len(voi)-1)):
            states[:,i+1] = sol.y[:,-1]

    #print(init_states)
    #print(sol.y)
    #print('size y=',np.shape(sol.y))
    #print(states)
    #print('size states=',np.shape(states))
    algebraic = computeAlgebraic(constants, states, voi)
    return (voi, states, algebraic)




def plot_model(voi, states, algebraic):
    """Plot variables against variable of integration"""
    import pylab
    (legend_states, legend_algebraic, legend_voi, legend_constants) = createLegends()
    pylab.figure(1)
    pylab.plot(voi,vstack((states,algebraic)).T)
    pylab.xlabel(legend_voi)
    #pylab.legend(legend_states + legend_algebraic, loc='best')
    pylab.show()

def plot_results(voi,states,algebraic):
    "This function is adopted by Hossein.Sharifi"
    (legend_states, legend_algebraic, legend_voi, legend_constants) = createLegends()

    from matplotlib import pyplot as plt
    import matplotlib.gridspec as gridspec
    import numpy as np


    #Ca legends
    #Ca_indicies_for_states=np.array([4,5,6,32,33,35,36,38])
    #Ca_indicies_for_states=np.array([32])
    Ca_indicies_for_states=np.array([32])
    #Ca_indicies_for_algebraic=np.array([46,47])

    size_Ca_for_states=len(Ca_indicies_for_states)
    #size_Ca_for_algebraic=len(Ca_indicies_for_algebraic)

    Ca_legend_states=[""] *size_Ca_for_states
    #Ca_legend_algebraic=[""] *size_Ca_for_algebraic

    Ca_states=array([[0.0] * len(voi)] * size_Ca_for_states)
    #Ca_algebraic=array([[0.0] * len(voi)] * size_Ca_for_algebraic)

    for i in range (0,size_Ca_for_states):
        Ca_legend_states[i]=legend_states[Ca_indicies_for_states[i]]
        Ca_states[i]=states[Ca_indicies_for_states[i]]
    #for i in range (0,size_Ca_for_algebraic):
    #    Ca_legend_algebraic[i]=legend_algebraic[Ca_indicies_for_algebraic[i]]
    #    Ca_algebraic[i]=algebraic[Ca_indicies_for_algebraic[i]]

    #print('size Ca_states=',np.shape(Ca_states))
    #print('size Ca_algebraic=',np.shape(Ca_algebraic))

    #states
    f=plt.figure(1,constrained_layout=True)
    f.set_size_inches([15,6])
    y_axis_states=vstack(Ca_states).T
    plt.plot(voi, y_axis_states)
    plt.xticks(np.arange(min(voi), max(voi)+1, 100.0))

    plt.xlabel(legend_voi)
    plt.ylabel('Ca_states (milimolar)')
    plt.legend(Ca_legend_states, bbox_to_anchor=(1.05, 1), \
                loc='best', borderaxespad=0.,fontsize='small')

    print("Saving Ca_states figure to")
    save_figure_to_file(f,"Shannon_Ca_states", dpi=None)

    #algebraic
    #f=plt.figure(2,constrained_layout=True)
    #y_axis_algebraic=vstack(Ca_algebraic).T
    #plt.plot(voi, vstack(y_axis_algebraic))
    #plt.xlabel(legend_voi)
    #plt.ylabel('algebraic')
    #plt.legend(Ca_legend_algebraic)
    #print("Saving Ca_algebraic figure to")
    #save_figure_to_file(f, "Ca_algebraic", dpi=None)
    #plt.show()

def save_figure_to_file(f,fname,dpi=None):
    "This function is adopted by Hossein.Sharifi"
    import os
    from skimage.io import imsave

    cwd=os.getcwd()
    filename=cwd + "/"+fname+".png"
    f.savefig(filename, dpi=dpi)

if __name__ == "__main__":
    #(voi, states, algebraic) = solve_model()
    (voi, states, algebraic) = solve_system()

    plot_results(voi, states, algebraic)
