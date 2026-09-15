#--------------------------------------------------------------------------------------------------------
# Solution of point defect balance equations in 1D
# Notes : 1- Equations are non-dimensionalized
#         2- Sinks(voids, bubbles, dislocations, grain boundaries) are uniformly distributed
#         3- Zero concentration is defined at domain(grain) boundaries
#         4- Defect generation is constant and uniform along boundary (Neutron Case)
#         5- UO2 parameters are used
#         6- Excess defect concentrations(C-Ceq) are calculated
#
# Reference textbook : https://doi.org/10.1007/978-1-4939-3438-6
# (Fundementals of Radiation Materials Science, Gary S. Was)
#--------------------------------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------------------------------
# Constant Parameters
# kB = 8.617e-5          #[ev/K] Boltzman Constant
# NA = 6.02214e23        #[1/mol] Avagadro's Number
# Rgas = ${fparse kB*NA} #[eV/mol-K] Universal Gas Constant : (1.987 cal/mol-K), (8.31446261815324 J/mol-K)
# pi = 3.14159265358979323846264338327950288

# [./ModelParameters]
T = 1000.0			# [Kelvin] Temperature [Kelvin]
d = 1e-06
Cv_ss = 2.3630683429819486e+29			# [1/m^3] Steady-State vacancy fraction [1/nm^3]
Ci_ss = 3854668552340487.5			# [1/m^3] Steady-State interstitial fraction [1/nm^3]
Cv_eq = 1140024855725.6003			# [1/m^3] Equilibrium vacancy fraction [1/nm^3]
Ci_eq = 0.22851719912638513			# [1/m^3] Equilibrium interstitial fraction [1/nm^3]
Cs = 1e+23
Di = 8.318668211156636e-16			# [m^2/s] Diffusivity for interstials [nm^2/s]
Dv = 1.6035945039527338e-19			# [m^2/s] Diffusivity for vacancies [nm^2/s]
K = 1e+24
Kiv = 5.489170003912054e-23			# [m^3/s] Recombination rate [nm3/s]
Kis = 5.488112057362333e-23			# [m^3/s] Interstitial sink reaction rate [nm^3/s]
Kvs = 1.0579465497205242e-26			# [m^3/s] Vacancy sink reaction rate [nm^3/s]
ki = 65973445725.38566			# [1/m^2] Interstitial sink strength [1/nm^2]
kv = 65973445725.38566			# [1/m^2] Vacancy sink strength [1/nm^2]
Va = 4.094533702e-29			# [m^3] Atomic volume [nm^3]
zi = 1.02             # Interstitial Sink bias zi~1.02 (2%)
# bias = 5.0            # Percent Production bias [%]
# eps = 5.0             # Percent FP production efficiency [%]
epsi = ${fparse 1-0.01*(27+0.01363636364*(T-700))}   #Interstitial production efficiency [fraction]
epsv = ${fparse 1-0.01*(39-0.01363636364*(T-700))}   #Vacancy production efficiency [fraction]
bias = ${fparse 100*(epsv/epsi-1)}  #Production bias factor [%]
tau1=${fparse 1/(K*Kiv)^0.5}   #Time Constant for Quasi Steady State
tau2=${fparse 1/(Kis*Cs)}      #Time Constant for interstitals to reach at sinks
tau3=${fparse 1/(Kvs*Cs)}      #Time Constant for vacancies to reach at sinks
tau4=${fparse Kis*Cs/(Kiv*K)}  #Time Constant for the transition between interstitials regimes (recombination dominates annihilation after this time)
# [../]

# [./Nondimensionalization Parameters]
nondim = 1	#  Control parameter for non-dimensionalization
ls = 1e-09	#  length scale
ts = ${fparse (1-nondim)+nondim*ls*ls/Di}   # time scale
omega= ${fparse (1-nondim)+nondim*Va}       # volume scale
# [../]

# [./Simulation Parameters]
l_cf = 1	#  unit conversion factor
sim_time = 1e20  # Simulation time
#h = ${ls}        # Element size
filename = 1D_UO2_pointdefects_ND# [../]

#--------------------------------------------------------------------------------------------------------
#  Modification is not necessary for rest of the file.
#--------------------------------------------------------------------------------------------------------
#
#---------------------------------------------------------------#
#                            Mesh                               #
#---------------------------------------------------------------#
[Mesh]
  type = GeneratedMesh
  dim = 1
  xmax = ${fparse d/ls}
  nx = 10000#${fparse int(d/h)}
[]
#---------------------------------------------------------------#
#                         Variables                             #
#---------------------------------------------------------------#
[Variables]
  [./xi]
  [../]
  [./xv]
  [../]
[]

[AuxVariables]
  [./ci]
  [../]
  [./cv]
  [../]
  [./v_supersaturation]
  [../]
  [./void_nucleation_rate]
  [../]
  # [./cv_analytical]
  #   type = MooseVariableFVReal
  # [../]
[]

#---------------------------------------------------------------#
#                          Kernels                              #
#---------------------------------------------------------------#

[Kernels]
  [./defect_generation_i]
    type = MaskedBodyForce
    variable = xi
    mask = source_i
  [../]
  [./defect_generation_v]
    type = MaskedBodyForce
    variable = xv
    mask = source_v
  [../]
  [./recombination_i]
    type = MatReaction
    variable = xi
    args = 'xv'
    mob_name = reaction_i
  [../]
  [./recombination_v]
    type = MatReaction
    variable = xv
    args = 'xi'
    mob_name = reaction_v
  [../]
  [./sink_reaction_i]
    type = MatReaction
    variable = xi
    mob_name = sink_i
  [../]
  [./sink_reaction_v]
    type = MatReaction
    variable = xv
    mob_name = sink_v
  [../]
  [./xi_diff]
    type = MatDiffusion
    variable = xi
    diffusivity = diff_i
  [../]
  [./xv_diff]
    type = MatDiffusion
    variable = xv
    diffusivity = diff_v
  [../]
  [./xi_time]
    type = TimeDerivative
    variable = xi
  [../]
  [./xv_time]
    type = TimeDerivative
    variable = xv
  [../]
[]

[AuxKernels]
  [./ci]
    type = ParsedAux
    variable = ci
    args = xi
    constant_names = 'omega l_cf'
    constant_expressions = '${omega} ${l_cf}'
    function = (l_cf^3)*xi/omega
  [../]
  [./cv]
    type = ParsedAux
    variable = cv
    args = xv
    constant_names = 'omega l_cf'
    constant_expressions = '${omega} ${l_cf}'
    function = (l_cf^3)*xv/omega
  [../]
  [./v_supersaturation]
    type = ParsedAux
    variable = v_supersaturation
    args = 'ci cv'
    constant_names = 'Di Dv Cv_eq'
    constant_expressions = '${Di} ${Dv} ${Cv_eq}'
    function = (Dv*cv-Di*ci)/(Dv*Cv_eq)
  [../]
  [./void_nucleation_rate]
    type = ParsedAux
    variable = void_nucleation_rate
    args = 'v_supersaturation'
    function = 'pow(v_supersaturation,5.41547)*exp(-14.6586)'
  [../]
  # [./cv_analytical]
  #   type = MaterialRealAux
  #   variable = cv_analytical
  #   property = mat_cv_analytical
  # [../]
[]

#---------------------------------------------------------------#
#                            BCs                                #
#---------------------------------------------------------------#
[BCs]
  [./xi_bc]
    type = DirichletBC
    variable = xi
    value = 0 #${fparse omega*Ci_eq}
    boundary = 'left right'
  [../]
  [./xv_bc]
    type = DirichletBC
    variable = xv
    value = 0 #${fparse omega*Cv_eq}
    boundary = 'left right'
  [../]
[]
#---------------------------------------------------------------#
#                            ICs                                #
#---------------------------------------------------------------#
[ICs]
  # [./xv] 
  #   type = RandomIC
  #   min = ${fparse 0.9*omega*Cv_eq}  # Equilibrium vacancy concentration
  #   max = ${fparse omega*Cv_eq}
  #   variable = xv
  # [../]
  # [./xi]
  #   type = RandomIC
  #   min = ${fparse 0.9*omega*Ci_eq}  # Equilibrium interstitial concentration
  #   max = ${fparse omega*Ci_eq}
  #   variable = xi
  # [../]
[]
#---------------------------------------------------------------#
#                         Materials                             #
#---------------------------------------------------------------#
[Materials]
  # Diffusion coefficients
  [./diff_i]
    type = ParsedMaterial
    f_name = diff_i
    material_property_names = 'Di w l'
    function = w*Di/l^2
  [../]
  [./diff_v]
    type = ParsedMaterial
    f_name = diff_v
    material_property_names = 'Dv w l'
    function = w*Dv/l^2
  [../]
  [./source_i]
    type = ParsedMaterial
    f_name = source_i
    material_property_names = 'K w omega epsi'
    function = epsi*w*omega*K
  [../]
  [./source_v]
    type = ParsedMaterial
    f_name = source_v
    material_property_names = 'K w omega bias epsv'
    function = epsv*w*omega*K
  [../]
  [./reaction_i]
    type = DerivativeParsedMaterial
    f_name = reaction_i
    args = 'xv'
    material_property_names = 'Kiv w omega'
    function = -w*Kiv*xv/omega
  [../]
  [./reaction_v]
    type = DerivativeParsedMaterial
    f_name = reaction_v
    args = 'xi'
    material_property_names = 'Kiv w omega'
    function = -w*Kiv*xi/omega
  [../]
  [./sink_i]
    type = DerivativeParsedMaterial
    f_name = sink_i
    material_property_names = 'Kis w omega Cs zi'
    function = -w*Kis*Cs*zi
  [../]
  [./sink_v]
    type = DerivativeParsedMaterial
    f_name = sink_v
    material_property_names = 'Kvs w omega Cs'
    function = -w*Kvs*Cs
  [../]
  [./dt]
    type = TimeStepMaterial
    prop_time = timestep
  [../]
  [./Constants]
    type = GenericConstantMaterial
    prop_names = 'w l omega Cs Di Dv K Kiv Kis Kvs ki kv d Va zi bias epsi epsv'
    prop_values = '${ts} ${ls} ${omega} ${Cs} ${Di} ${Dv} ${K} ${Kiv} ${Kis} ${Kvs} ${ki} ${kv} ${d} ${Va} ${zi} ${bias} ${epsi} ${epsv}'
  [../]
  # [./cv_analytical]
  #   type = ParsedAux
  #   f_name = mat_cv_analytical
  #   material_property_names = 'source_v diff_v sink_v'
  #   function = (source_v/sink_v)*(1-exp(-sqrt(x*sink_v/diff_v)))
  # [../]
[]
#---------------------------------------------------------------#
#                       Post Processors                         #
#---------------------------------------------------------------#
[Postprocessors]
  [./xv]
    type = ElementAverageValue
    variable = xv
  [../]
  [./xi]
    type = ElementAverageValue
    variable = xi
  [../]
  [./cv]
    type = ParsedPostprocessor
    pp_names = 'xv omega'
    constant_names = 'l_cf'
    constant_expressions ='${l_cf}'
    function = (l_cf^3)*xv/omega
    # outputs = 'exodus csv'
  [../]
  [./ci]
    type = ParsedPostprocessor
    pp_names = 'xi omega'
    constant_names = 'l_cf'
    constant_expressions = '${l_cf}'
    function = (l_cf^3)*xi/omega
    # outputs = 'exodus csv'
  [../]
  [ci_abs]
    type = ParsedPostprocessor
    pp_names = 'ci ci_eq'
    function = ci+ci_eq
    outputs = 'csv'
  []
  [cv_abs]
    type = ParsedPostprocessor
    pp_names = 'cv cv_eq'
    function = cv+cv_eq
    outputs = 'csv'
  []
  [./total_xv]
    type = ElementIntegralVariablePostprocessor
    variable = xv
    outputs = 'csv'
  [../]
  [./total_xi]
    type = ElementIntegralVariablePostprocessor
    variable = xi
    outputs = 'csv'
  [../]
  [./right_jvx]
    type = SideDiffusiveFluxAverage
    variable = xv
    boundary = right
    diffusivity = diff_v
    # outputs = 'exodus csv'
  [../]
  [./left_jvx]
    type = SideDiffusiveFluxAverage
    variable = xv
    boundary = left
    diffusivity = diff_v
    outputs = 'csv'
  [../]
  [./right_jix]
    type = SideDiffusiveFluxAverage
    variable = xi
    boundary = right
    diffusivity = diff_i
    # outputs = 'exodus csv'
  [../]
  [./left_jix]
    type = SideDiffusiveFluxAverage
    variable = xi
    boundary = left
    diffusivity = diff_i
    outputs = 'csv'
  [../]
  [./t]
    type = ParsedPostprocessor
    pp_names = 'w'
    use_t = true
    function = t*w
    outputs = 'csv'
  [../]
  [./t_hours]
    type = ParsedPostprocessor
    pp_names = 'w'
    use_t = true
    function = t*w/60/60
    outputs = 'csv'
  [../]
  [./t_days]
    type = ParsedPostprocessor
    pp_names = 'w'
    use_t = true
    function = t*w/60/60/24
    outputs = 'csv'
  [../]
  [./w]
    type = ElementAverageMaterialProperty
    mat_prop = w
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  [../]
  [./l]
    type = ElementAverageMaterialProperty
    mat_prop = l
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  [../]
  [./omega]
    type = ElementAverageMaterialProperty
    mat_prop = omega
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  [../]
  [./d]
    type = ElementAverageMaterialProperty
    mat_prop = d
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  [../]
  [bias]
    type = ParsedPostprocessor
    pp_names = 'source_i source_v'
    constant_names = bias
    constant_expressions = ${bias}
    function = bias
    outputs = 'exodus csv'
  []
  [./source_i]
    type = ElementAverageMaterialProperty
    mat_prop = source_i
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  [../]
  [./source_v]
    type = ElementAverageMaterialProperty
    mat_prop = source_v
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  [../]
  [./diff_i]
    type = ElementAverageMaterialProperty
    mat_prop = diff_i
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  [../]
  [./diff_v]
    type = ElementAverageMaterialProperty
    mat_prop = diff_v
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  [../]
  [./reaction_i]
    type = ElementAverageMaterialProperty
    mat_prop = reaction_i
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  [../]
  [./reaction_v]
    type = ElementAverageMaterialProperty
    mat_prop = reaction_v
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  [../]
  [./sink_i]
    type = ElementAverageMaterialProperty
    mat_prop = sink_i
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  [../]
  [./sink_v]
    type = ElementAverageMaterialProperty
    mat_prop = sink_v
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  [../]
  [./v_supersaturation]
    type = ElementAverageValue
    variable = v_supersaturation
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  [../]
  [./void_nucleation_rate]
    type = ElementAverageValue
    variable = void_nucleation_rate
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  [../]
  [./kgbi]
    type = ParsedPostprocessor
    pp_names = 'ci cv diff_i diff_v omega w l left_jix left_jvx d'
    function = (6*left_jix/(omega*w/l))/(d*ci*diff_i*l^2/w)
    execute_on = 'TIMESTEP_END'
    outputs = 'csv'
  [../]
  [./kgbv]
    type = ParsedPostprocessor
    pp_names = 'ci cv diff_i diff_v omega w l left_jix left_jvx d'
    function = (6*left_jvx/(omega*w/l))/(d*cv*diff_v*l^2/w)
    execute_on = 'TIMESTEP_END'
    outputs = 'csv'
  [../]
  [./Dratio]
    type = ParsedPostprocessor
    pp_names = 'diff_i diff_v'
    function = diff_i/diff_v
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = 'csv'
  [../]
  [ci_ss]
    type = ParsedPostprocessor
    pp_names = 'ci'
    constant_names = ci_ss
    constant_expressions = ${Ci_ss}
    function = ci_ss
    outputs = 'csv'
  []
  [cv_ss]
    type = ParsedPostprocessor
    pp_names = 'cv'
    constant_names = cv_ss
    constant_expressions = ${Cv_ss}
    function = cv_ss
    outputs = 'csv'
  []
  [ci_eq]
    type = ParsedPostprocessor
    pp_names = 'ci'
    constant_names = ci_eq
    constant_expressions = ${Ci_eq}
    function = ci_eq
    outputs = 'csv'
  []
  [cv_eq]
    type = ParsedPostprocessor
    pp_names = 'cv'
    constant_names = cv_eq
    constant_expressions = ${Cv_eq}
    function = cv_eq
    outputs = 'csv'
  []
  [tau1]
    type = ParsedPostprocessor
    pp_names = 'ci'
    constant_names = tau1
    constant_expressions = ${tau1}
    function = tau1
    outputs = 'csv'
  []
  [tau2]
    type = ParsedPostprocessor
    pp_names = 'ci'
    constant_names = tau2
    constant_expressions = ${tau2}
    function = tau2
    outputs = 'csv'
  []
  [tau3]
    type = ParsedPostprocessor
    pp_names = 'ci'
    constant_names = tau3
    constant_expressions = ${tau3}
    function = tau3
    outputs = 'csv'
  []
  [tau4]
    type = ParsedPostprocessor
    pp_names = 'ci'
    constant_names = tau4
    constant_expressions = ${tau4}
    function = tau4
    outputs = 'csv'
  []
[]
[VectorPostprocessors]
  [./x_direc]
   type =  LineValueSampler
    start_point = '0 0 0'
    end_point = '${fparse d/ls} 0 0'
    variable = 'xi xv ci cv v_supersaturation void_nucleation_rate'
    num_points = 1001
    sort_by =  id
  [../]
[]
#---------------------------------------------------------------#
#                       Preconditioning                         #
#---------------------------------------------------------------#
[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]
#---------------------------------------------------------------#
#                        Executioner                            #
#---------------------------------------------------------------#
[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu  superlu_dist'
  resid_vs_jac_scaling_param = 0.5
  scheme = bdf2
  automatic_scaling = true
  line_search = none
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10  # 1/100 of lowest  equilibrium defect concentration
  l_max_its = 50
  l_tol = 1e-4
  dt = 1
  start_time = 0
  end_time = ${sim_time}  # 1 day = 86400 sec
  dtmax = 1e15
  dtmin = 1e-3
  steady_state_detection = true
  steady_state_tolerance = 1e-15
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 7
    growth_factor = 1.8
    cutback_factor = 0.1
    dt = 1e-3
  []
[]
#---------------------------------------------------------------#
#                          Outputs                              #
#---------------------------------------------------------------#
[Outputs]
  file_base = ${filename}
  [./exodus]
    type = Exodus
    enable = false #exodus
    output_material_properties = true
    execute_postprocessors_on = 'INITIAL TIMESTEP_END FAILED'
    # interval = 100 #exodus
  [../]
  [./csv]
    type = CSV
    execute_postprocessors_on = 'INITIAL TIMESTEP_END FAILED'
    execute_vector_postprocessors_on = 'FINAL FAILED'
  [../]
[]
#--------------------------------------------------------------------------------------------------------