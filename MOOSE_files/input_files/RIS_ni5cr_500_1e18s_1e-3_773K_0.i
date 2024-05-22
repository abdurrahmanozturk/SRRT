# Ni-Cr in 1D
# l_scale = 1e-9 m
# Equations are non-dimensionalized

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 500
  #ny = 64
  #nz = 0
  xmax = 500
  #ymax = 256
[]

[Variables]
  [./Xv]
  [../]
  [./Xi]
  [../]
  [./X_cr] #Cr
  [../]
[]

[AuxVariables]
  [./X_ni]
  [../]
  #T = * 773K values
  [./dcri] #chromium
    initial_condition = 1.2380657478440465e-09
  [../]
  [./dfei] #nickel
    initial_condition = 5.577318586890966e-10
  [../]
  [./dcrv] #chromium
    initial_condition = 1.5208619674807823e-11
  [../]
  [./dfev] #nickel
    initial_condition = 5.057688432327799e-12
  [../]
  [./Xs] #Uniform sink concentration
    initial_condition = 1.206e-11
  [../]
[]

[ICs]
  [./Xv] # Concentration of vacancy
    type = RandomIC
    variable = Xv
    min = 2.0e-12
    max = 2.2e-12
    seed = 11
[../]

[./Xi] # Concentration of interstials
  type = RandomIC
  variable = Xi
  min = 8.2e-27
  max = 8.4e-27
  seed = 11
[../]

[./X_cr] # Concentration of Cr
  type = ConstantIC
   value = 0.05
  variable = X_cr
  #min = 0.08
  #max = 0.1
  #seed = 11
[../]

[]

[BCs]

  [./Xv]
  type = DirichletBC
  preset = false
  variable = 'Xv'
  boundary = 'left'
  value = 2.133878404598678e-12
  [../]

  [./Xi]
  type = DirichletBC
  preset = false
  variable = 'Xi'
  boundary = 'left'
  value = 8.316339677291549e-27
  [../]
[]

[AuxKernels]

  [./X_ni]
  type = ParsedAux
  variable = X_ni
  args = 'X_cr'
  function = '1-X_cr'
  [../]
[]

[Kernels]

                            #### Xv_equation ###

  [./Xv_dot]
   type = TimeDerivative
    variable = 'Xv'
  [../]

  [./Xv_dcrv_dfev]
    type = MatDiffusion
    variable = Xv
    D_name = dcrv_dfev
    v = X_cr
    args = Xv
  [../]

  [./Xv_Dv]
    type = MatDiffusion
    variable = Xv
    D_name = Xv_Dv
    args = 'X_cr'
  [../]

  [./Xv_source_K0]
      type = MaskedBodyForce
      variable = Xv
      mask = Xv_source_K0
      args = ' '
  [../]

  [./Xv_KivXiXv]
      type = MatReaction
      variable = Xv
      mob_name = Xv_KivXiXv
      args = 'Xi X_cr'
  [../]

  [./Sink_v]
    type = MatReaction
    variable = Xv
    args = 'Xs X_cr'
    mob_name = sink_v
  [../]

                          ###   Xi_equation  ####

  [./Xi_dot]
    type = TimeDerivative
    variable = 'Xi'
  [../]

  [./Xi_dcri_dfei]
    type = MatDiffusion
    variable = Xi
    D_name = dcri_dfei
    v = X_cr
    args = Xi
  [../]

  [./Xi_Di]
    type = MatDiffusion
    variable = Xi
    D_name = Xi_Di
    args = 'X_cr'
  [../]

  [./Xi_source_K0]
    type = MaskedBodyForce
    variable = Xi
    mask = Xi_source_K0
    args = ' '
  [../]

  [./Xi_KivXiXv]
        type = MatReaction
        variable = Xi
        mob_name = Xi_KivXiXv
        args = 'Xv X_cr'
  [../]

  [./Sink_i]
      type = MatReaction
      variable = Xi
      args = 'Xs X_cr'
      mob_name = sink_i
  [../]

                        #### X_Cr_equation  ####

  [./X_cr_dot]
    type = TimeDerivative
    variable = 'X_cr'
  [../]

  [./X_cr_Dcr]
    type = MatDiffusion
    variable = X_cr
    D_name = X_cr_Dcr
    args = 'Xv Xi'
  [../]

  [./X_cr_dcri]
    type = MatDiffusion
    variable = X_cr
    D_name = dcri_bar
    v = Xi
    args = X_cr
  [../]

  [./X_cr_dcrv]
    type = MatDiffusion
    variable = X_cr
    D_name = dcrv_bar
    v = Xv
    args = X_cr
  [../]

[]

[Materials]

                          ##### Xv_equation_material###

        [./Xv_dcrv_dfev]
          type = DerivativeParsedMaterial
          constant_names = 'chi' #chi: thermodynamic_factor (unitless)
          constant_expressions = '1'
          f_name = dcrv_dfev
          function = -(chi)*Xv*(dcrv-dfev)/(dcri)
          args = 'Xv dcrv dcri dfev'
        [../]

        [./Xv_Dv]
          type = DerivativeParsedMaterial
          f_name = Xv_Dv
          function = (((dcrv-dfev)*X_cr+dfev)/dcri)
          args = 'X_cr dcrv dcri dfev'
        [../]

        [./Xv_source_K0]
          type = DerivativeParsedMaterial
          constant_names = 'K0 l_scale'
          constant_expressions = '1e-3 1e-9'
          #production bias
          f_name = Xv_source_K0
          function = K0*((l_scale*l_scale)/dcri)
          args = 'dcri'
        [../]

        [./Xv_KivXiXv]
          type = DerivativeParsedMaterial
          constant_names = 'l_scale omega'
          constant_expressions = '1e-9 1.206e-29'
          f_name = Xv_KivXiXv
          function = -4*3.14*10*(0.352e-9)*(((dcri-dfei)*X_cr+dfei)+((dcrv-dfev)*X_cr+dfev))*((l_scale*l_scale)/(dcri*omega))*Xi
          args = 'Xi X_cr dcri dfei dfev dcrv'
        [../]

        [./Sink_v]
          type = DerivativeParsedMaterial
          constant_names = 'l_scale omega '
          constant_expressions = ' 1e-9 1.206e-29'
          f_name = sink_v
          args = 'Xs X_cr dcri dfei dcrv dfev'
          function = -4*3.14*10*(0.352e-9)*((dcrv-dfev)*X_cr+dfev)*((l_scale*l_scale)/(dcri*omega))*Xs
        [../]

                          ###  Xi_equation_material #####

        [./Xi_dcri_dfei]
          type = DerivativeParsedMaterial
          constant_names = 'chi omega'
          constant_expressions = '1 1.206e-29'
          f_name = dcri_dfei
          function = chi*Xi*(dcri-dfei)/dcri
          args = 'Xi dcri dfei'
        [../]

        [./Xi_Di]
          type = DerivativeParsedMaterial
          f_name = Xi_Di
          function = (((dcri-dfei)*X_cr+dfei)/dcri)
          args = 'X_cr dcri dfei'
        [../]

        [./Xi_source_K0]
          type = DerivativeParsedMaterial
          constant_names = 'K0 l_scale '
          constant_expressions = '1.0e-3 1e-9'
          f_name = Xi_source_K0
          function = K0*((l_scale*l_scale)/dcri)
          args = 'dcri'
        [../]

        [./Xi_KivXiXv]
          type = DerivativeParsedMaterial
          constant_names = 'l_scale omega '
          constant_expressions = ' 1e-9 1.206e-29'
          f_name = Xi_KivXiXv
          function = -4*3.14*10*(0.352e-9)*(((dcri-dfei)*X_cr+dfei)+((dcrv-dfev)*X_cr+dfev))*((l_scale*l_scale)/(dcri*omega))*Xv
          args = 'Xv X_cr dcri dfei dcrv dfev'
        [../]

        [./sink_i]
          type = DerivativeParsedMaterial
          constant_names = 'l_scale omega '
          constant_expressions = ' 1e-9 1.206e-29'
          f_name = sink_i
          args = 'Xs X_cr dcri dfei dcrv dfev'
          function = -4*3.14*10*(0.352e-9)*((dcri-dfei)*X_cr+dfei)*((l_scale*l_scale)/(dcri*omega))*Xs
        [../]


                            #### X_Cr_equation_material  ####
       [./X_cr_Dcr]
          type = DerivativeParsedMaterial
          constant_names = 'chi'
          constant_expressions = '1'
          f_name = X_cr_Dcr
          function = chi*(dcrv*Xv+dcri*Xi)/dcri
          args = 'dcri dcrv Xi Xv'
       [../]

       [./X_cr_dcri]
         type = DerivativeParsedMaterial
         f_name = dcri_bar
         function = X_cr*dcri/dcri
         args = 'X_cr dcri'
       [../]

       [./X_cr_dcrv]
         type = DerivativeParsedMaterial
         f_name = dcrv_bar
         function = -X_cr*dcrv/dcri
         args = 'X_cr dcri dcrv'
       [../]


[]

[Postprocessors]

        [./tot_Xv]
          type = ElementIntegralVariablePostprocessor
          variable = Xv
        [../]
        [./tot_Xi]
          type = ElementIntegralVariablePostprocessor
          variable = Xi
        [../]
        [./tot_X_cr]
          type = ElementIntegralVariablePostprocessor
          variable = X_cr
        [../]
        [./average_Xi]
          type = ElementAverageValue
          variable = Xi
        [../]
        [./average_Xv]
          type = ElementAverageValue
          variable = Xv
        [../]
        [./average_X_cr]
          type = ElementAverageValue
          variable = X_cr
        [../]
        [./left_X_cr]
          type = PointValue
          point = '0.0 0.0 0.0'
          variable = X_cr
        [../]
        [./right_X_cr]
          type = PointValue
          point = '500 0.0 0.0'
          variable = X_cr
        [../]
[]

[VectorPostprocessors]
        [./x_direc]
         type =  LineValueSampler
          start_point = '0 0 0'
          end_point = '500 0 0'
          variable = 'X_cr Xv Xi X_ni'
          num_points = 501
          sort_by =  id
        [../]
[]

[Preconditioning]
        [./SMP]
          type = SMP
          full = true
        [../]
[]

[Executioner]
        # Preconditioned JFNK (default)
        scheme = BDF2
        type = Transient
        nl_max_its = 10
        solve_type = NEWTON
      #  petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_tol -pc_hypre_boomeramg_max_iter'
      #  petsc_options_value = 'hypre boomeramg 1e-4 20'
         petsc_options_iname = '-pc_type'
         petsc_options_value = 'lu'
        l_max_its = 15 #max linear iterations default 10000
        l_tol = 1.0e-3 #linear tolerance
        nl_rel_tol = 1.0e-6 #nonlinear relative tolerance default 1e-8
        start_time = 0.0
        nl_abs_tol = 1e-14
        num_steps = 150000
        #automatic_scaling =  true
        steady_state_detection = true
         steady_state_tolerance = 1e-16
        [./TimeStepper]
          type = IterationAdaptiveDT
          cutback_factor = .75
          dt = 1.0
          growth_factor = 1.2
          optimal_iterations = 7
        [../]
[]
[Outputs]
        csv = true
        exodus = true
        interval = 1
[]
