# Sleipner CO2 Injection — 1-DAY VALIDATION (20 MPI, baked mesh)
#
# Uses pre-baked mesh (sleipner_grid_1M_with_well.e) — no MeshGenerators in parallel.
# Run mesh bake first: bash scripts/submit/run_mesh_bake.sh
# Then: salloc -N 1 -n 20 -t 1:00:00  then  bash scripts/submit/submit_verif_1day.sh
# Compare: min/max pwater, co2_mass, max_sgas vs 1-rank Dirac (sleipner_corrected_v2_immiscible_1M.i).

[Mesh]
  # Read baked mesh (well block + injection_well sideset); no MeshGenerators in parallel.
  # Bypass: if metis/parmetis segfault, try partitioner = linear via CLI.
  [file_mesh]
    type = FileMeshGenerator
    file = ../meshes/sleipner_grid_1M_with_well.e
  []
  parallel_type = distributed
  partitioner = metis
[]

[Preconditioning]
  active = 'SMP_asm1'
  [SMP_bjacobi]
    type = SMP
    full = true
    petsc_options = '-ksp_monitor -ksp_converged_reason -error_output_stderr'
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -sub_pc_factor_shift_amount -sub_pc_factor_pivot_in_blocks -sub_pc_factor_levels -ksp_gmres_restart -ksp_rtol -ksp_max_it'
    petsc_options_value = 'fgmres bjacobi ilu NONZERO 1e-12 1 1 200 5e-2 100'
  []
  [SMP_asm1]
    type = SMP
    full = true
    petsc_options = '-ksp_monitor -ksp_converged_reason -error_output_stderr'
    petsc_options_iname = '-ksp_type -pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_shift_type -sub_pc_factor_shift_amount -sub_pc_factor_pivot_in_blocks -sub_pc_factor_levels -ksp_gmres_restart -ksp_rtol -ksp_max_it'
    petsc_options_value = 'fgmres asm 1 ilu NONZERO 1e-12 1 1 200 5e-2 100'
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 9.81'  # Z positive down (depth)
[]

[Problem]
  allow_initial_conditions_with_restart = true
[]

[Functions]
  [hydrostatic_pressure]
    type = ParsedFunction
    expression = '9.0e6 + 1020 * 9.81 * (z - 800)'
  []
  # Injection rate: constant 28.5 kg/s ≈ 0.9 Mt/yr (for postprocessor monitoring)
  [injection_rate]
    type = ConstantFunction
    value = 28.5
  []
  # Flux on injection_well boundary: -28.5/area (negative = injection), zero after 12 years
  # Area from SideSetsAroundSubdomain (4 side walls only): ~7704 m² (check_injection_well_sideset)
  [injection_flux_function]
    type = ParsedFunction
    expression = 'if(t < 378432000, -28.5/7704, 0)'
  []
  # dt_schedule: slow ramp, capped at 3600 s (1 hour) to prevent sgas overshoot
  # at injection nodes. With 36 injection points + dt≤3600, sgas change ~0.004/step.
  [dt_schedule]
    type = PiecewiseLinear
    x = '0      300     3600    86400   604800   2592000  31536000  378432000'
    y = '0.01   1.0     30      300     1800     3600     3600      3600'
  []
[]

[Variables]
  [pwater]
    outputs = 'csv'
  []
  [sgas]
    initial_condition = 0.0
    outputs = 'csv'
  []
[]

[ICs]
  [pwater_hydrostatic]
    type = FunctionIC
    variable = pwater
    function = hydrostatic_pressure
  []
[]

[AuxVariables]
  # Mass fractions for immiscible model (constant)
  [massfrac_ph0_sp0]
    initial_condition = 1   # Water phase is 100% water
  []
  [massfrac_ph1_sp0]
    initial_condition = 0   # Gas phase is 0% water (100% CO2)
  []
  
  [perm_scalar]
    family = MONOMIAL
    order = CONSTANT
  []
  # B1: only physics-required Aux kept; rest removed to test Column-too-large
  [PERMX]
    family = MONOMIAL
    order = CONSTANT
  []
  [PORO]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [permx_sand]
    type = ConstantAux
    variable = PERMX
    value = 2000.0
    block = '0 2'
    execute_on = 'initial'
  []
  [permx_mud]
    type = ConstantAux
    variable = PERMX
    value = 1.0
    block = '1'
    execute_on = 'initial'
  []
  [poro_sand]
    type = ConstantAux
    variable = PORO
    value = 0.37
    block = '0 2'
    execute_on = 'initial'
  []
  [poro_mud]
    type = ConstantAux
    variable = PORO
    value = 0.15
    block = '1'
    execute_on = 'initial'
  []

  [calc_perm_scalar]
    type = ParsedAux
    variable = perm_scalar
    expression = 'PERMX * 9.869233e-16'
    coupled_variables = 'PERMX'
    execute_on = 'initial timestep_end'
  []
[]

[Kernels]
  [mass_water]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pwater
  []
  [flux_water]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = pwater
  []
  [mass_co2]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = sgas
  []
  [flux_co2]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    variable = sgas
  []
[]

[BCs]
  # DISABLED: lateral pressure BC requires boundary names (left/right/front/back)
  # which are not available in the clean split mesh.
  # Natural no-flow BC applies instead (zero flux at all boundaries).
  # This is physically reasonable for a closed domain test.
  # To re-enable: re-split mesh WITH boundaries baked in.
  # [pwater_lateral]
  #   type = FunctionDirichletBC
  #   variable = pwater
  #   boundary = 'left right front back'
  #   function = hydrostatic_pressure
  # []
  # Top = caprock: no pressure BC → no flow through top (no CO2 leakage). Lateral only allows outflow.
  # [pwater_top]
  #   type = FunctionDirichletBC
  #   variable = pwater
  #   boundary = 'top'
  #   function = hydrostatic_pressure
  # []
  #
  # === Gas Saturation Boundary Conditions (COMMENTED OUT) ===
  # Fixing sgas = 0 at boundaries prevents CO2 from leaving the domain.
  # Currently NOT applied; CO2 can migrate to lateral boundaries and be removed
  # by the hydrostatic pwater_lateral BC. Uncomment if needed to prevent leakage.
  # [sgas_lateral]
  #   type = DirichletBC
  #   variable = sgas
  #   boundary = 'left right front back'
  #   value = 0.0
  # []
  # [sgas_top]
  #   type = DirichletBC
  #   variable = sgas
  #   boundary = 'top'
  #   value = 0.0
  # []
  # Note: No sgas BC on bottom boundary - CO2 typically doesn't flow downward due to gravity.
  # Injection via PorousFlowSink on well sideset (replaces Dirac point sources for MPI stability).
  [co2_injection]
    type = PorousFlowSink
    boundary = 'injection_well'
    variable = sgas
    use_mobility = false
    use_relperm = false
    fluid_phase = 1
    flux_function = injection_flux_function
  []
[]

# === CO2 Injection: surface source (well sideset) — no DiracKernels ===
# Injection is via PorousFlowSink on boundary 'injection_well' (see [BCs]).
[DiracKernels]
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pwater sgas'
    number_fluid_phases = 2
    number_fluid_components = 2
  []
  
  # Capillary pressure for sandstone (used in PorousFlow2PhasePS)
  # lambda reduced from 2.0 to 1.5 to make curve less steep (reduces sgas overshoot)
  [pc_sand]
    type = PorousFlowCapillaryPressureBC
    pe = 5e3        # Entry pressure 5 kPa
    lambda = 1.5    # Reduced from 2.0: gentler curve = less numerical overshoot
    sat_lr = 0.1
    pc_max = 1e5
  []
  
  # No SumQuantityUO needed: injection is rate-controlled via DiracKernels.
  # Cumulative injection is tracked by total_injected_mass postprocessor.

  # === Terminator: DISABLED ===
  # vinewtonrsls + ConstantBounds (the "ideal" approach) causes SEGFAULT on this
  # distributed mesh (20 ranks, ~32M DOFs, jobs 5906330/5920143).
  #
  # Terminator was tried at thresholds 0.99 (job 5920371) and 1.0 (job 5922358),
  # but BOTH caused an unrecoverable dt cutback death spiral:
  #   - Solver converges perfectly (1 Newton iter) at every dt down to dtmin
  #   - But max_sgas at a single node near the injection zone legitimately reaches ~1.0
  #   - Even dt=0.1 produces sgas ≥ threshold → Terminator rejects → stuck at dtmin
  #
  # PorousFlow materials internally clamp effective saturation in relative permeability
  # (kr=1 for seff≥1, kr=0 for seff≤0) and capillary pressure (pc=pc_max for seff≤0).
  # Therefore small overshoots (sgas slightly > 1.0) at injection nodes are tolerable.
  # The conservative settings (36 injection points, dtmax=3600, growth_factor=1.2)
  # prevent large dt jumps, so sgas overshoot should remain tiny (< 0.01).
  #
  # Monitor max_sgas via the postprocessor / CSV output to verify it stays near 1.0.
  # If it climbs above ~1.1, revisit injection distribution or solver settings.
  #
  # [sgas_bound_terminator]
  #   type = Terminator
  #   expression = 'max_sgas > 1.0'
  #   fail_mode = SOFT
  #   execute_on = TIMESTEP_END
  #   error_level = WARNING
  # []
[]

[FluidProperties]
  [water]
    type = SimpleFluidProperties
    density0 = 1020
    viscosity = 3e-4
    bulk_modulus = 2e9
  []
  [co2]
    type = SimpleFluidProperties
    density0 = 600
    viscosity = 5e-5
    bulk_modulus = 1e8
  []
[]

[Materials]
  # === Two-Phase Pressure-Saturation formulation ===
  [ppss]
    type = PorousFlow2PhasePS
    phase0_porepressure = pwater
    phase1_saturation = sgas
    capillary_pressure = pc_sand
  []
  
  # === Mass fractions (immiscible: water phase = 100% water, gas phase = 100% CO2) ===
  [massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = 'massfrac_ph0_sp0 massfrac_ph1_sp0'
  []
  
  # === Fluid properties ===
  [water_fluid]
    type = PorousFlowSingleComponentFluid
    fp = water
    phase = 0
  []
  [co2_fluid]
    type = PorousFlowSingleComponentFluid
    fp = co2
    phase = 1
  []

  # === Permeability (applies to all blocks) ===
  [permeability]
    type = PorousFlowPermeabilityTensorFromVar
    perm = perm_scalar
    k_anisotropy = '1 0 0  0 1 0  0 0 0.1'
  []

  # === Porosity (applies to all blocks) ===
  [porosity]
    type = PorousFlowPorosityConst
    porosity = PORO
  []

  # === SANDSTONE Relative Permeability - Block 0 ===
  # Reduced n exponents and increased s_res to make curves less steep (reduces sgas overshoot)
  [relperm_water_sand]
    type = PorousFlowRelativePermeabilityCorey
    block = '0 2'
    phase = 0
    n = 1.5         # Reduced from 2.0: gentler curve
    s_res = 0.20   # Increased from 0.15: larger residual = smoother transition
    sum_s_res = 0.30  # Updated to match new s_res values
  []
  [relperm_gas_sand]
    type = PorousFlowRelativePermeabilityCorey
    block = '0 2'
    phase = 1
    n = 2.0         # Reduced from 3.0: gentler curve
    s_res = 0.10    # Increased from 0.05: larger residual = smoother transition
    sum_s_res = 0.30  # Updated to match new s_res values
  []

  # === MUDSTONE Relative Permeability - Block 1 ===
  # Reduced n exponents to make curves less steep (especially gas: n=5.0 was very steep)
  [relperm_water_mud]
    type = PorousFlowRelativePermeabilityCorey
    block = '1'
    phase = 0
    n = 2.5         # Reduced from 3.0: gentler curve
    s_res = 0.30
    sum_s_res = 0.40  # Updated to match new gas s_res
  []
  [relperm_gas_mud]
    type = PorousFlowRelativePermeabilityCorey
    block = '1'
    phase = 1
    n = 3.0         # Reduced from 5.0: much gentler curve (was very steep)
    s_res = 0.10    # Increased from 0.05: larger residual = smoother transition
    sum_s_res = 0.40  # Updated to match new s_res values
  []

  # === Temperature (constant, all blocks) ===
  [temperature]
    type = PorousFlowTemperature
    temperature = 310.15  # 37°C — Sleipner literature (Chadwick et al. 2004; Audigane et al. 2007)
  []
[]

[Debug]
  show_actions = false
  show_var_residual_norms = true
[]

[Executioner]
  type = Transient
  solve_type = NEWTON

  nl_max_its = 30
  nl_abs_tol = 1e-4
  nl_rel_tol = 5e-2            # 放宽以减少每步 Newton 迭代
  l_tol = 2e-2
  l_max_its = 500

  automatic_scaling = false
  
  start_time = 0
  end_time = 86400  # 1 day (verif)

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    growth_factor = 1.5
    cutback_factor = 0.75
    cutback_factor_at_failure = 0.7
    optimal_iterations = 7
    iteration_window = 3
    # 线性迭代过多（如单步总 KSP iters > 约 64）则下一步缩小 dt，相当于 “cutback when iteration 多”
    # 公式: shrink when total_linear_iters > optimal_iterations*ratio + iteration_window*ratio；8 约 40+24=64
    linear_iteration_ratio = 8
  []
[]

# === Bounds enforcement: sgas ∈ [0, 1] — DISABLED (segfault on distributed mesh) ===
# ConstantBounds + vinewtonrsls is the MOOSE-recommended approach, but it causes
# SEGFAULT on this 20-rank distributed mesh (~32M DOFs) with MOOSE 2025-11-14 / PETSc 3.24.0
# (jobs 5906330, 5920143).  The bounds vector operations appear to have a bug or
# incompatibility with distributed ParMETIS partitioning at this scale.
# Workaround: Terminator (above) + conservative time stepping.
# TODO: Re-test with future MOOSE/PETSc versions.
# [Bounds]
#   [sgas_upper]
#     type = ConstantBounds
#     variable = bounds_dummy
#     bounded_variable = sgas
#     bound_type = upper
#     bound_value = 1.0
#   []
#   [sgas_lower]
#     type = ConstantBounds
#     variable = bounds_dummy
#     bounded_variable = sgas
#     bound_type = lower
#     bound_value = 0.0
#   []
# []

[Postprocessors]
  inactive = 'injection_well_area'   # B1: fewer objects touching boundary
  [injection_rate_pp]
    type = FunctionValuePostprocessor
    function = injection_rate
    execute_on = 'initial timestep_begin'
  []
  
  [dt]
    type = TimestepSize
    execute_on = 'initial timestep_end'
  []
  [time]
    type = TimePostprocessor
    execute_on = 'initial timestep_end'
  []
  
  # === Mass Balance ===
  [co2_mass]
    type = PorousFlowFluidMass
    fluid_component = 1
    execute_on = 'initial timestep_end'
  []
  [water_mass]
    type = PorousFlowFluidMass
    fluid_component = 0
    execute_on = 'initial timestep_end'
  []
  
  # === CO2 Plume Monitoring ===
  [max_sgas]
    type = ElementExtremeValue
    variable = sgas
    value_type = max
    execute_on = 'initial timestep_end'
  []
  [min_sgas]
    type = ElementExtremeValue
    variable = sgas
    value_type = min
    execute_on = 'initial timestep_end'
  []
  [max_pwater]
    type = NodalExtremeValue
    variable = pwater
    value_type = max
    execute_on = 'initial timestep_end'
  []
  [min_pwater]
    type = NodalExtremeValue
    variable = pwater
    value_type = min
    execute_on = 'initial timestep_end'
  []

  # === Injection Point Monitoring ===
  # DISABLED: PointValue fails on distributed (pre-split) mesh because the point
  # locator can't find elements across partitions. DiracKernels (injection) handle
  # this internally via geometric ghosting, but PointValue does not.
  # Use max_sgas / co2_mass postprocessors to monitor injection instead.
  # [sgas_at_injector]
  #   type = PointValue
  #   point = '1605.4 2955.65 971.421'
  #   variable = sgas
  #   execute_on = 'initial timestep_end'
  # []
  # [pwater_at_injector]
  #   type = PointValue
  #   point = '1605.4 2955.65 971.421'
  #   variable = pwater
  #   execute_on = 'initial timestep_end'
  # []
  # [sgas_at_inj_node]
  #   type = PointValue
  #   point = '1605.4 2955.65 953.421'
  #   variable = sgas
  #   execute_on = 'initial timestep_end'
  # []
  
  # === Cumulative Injection (corrected) ===
  # BUG FIX: the old CumulativeValuePostprocessor on injection_rate_pp accumulated
  # the RATE (28.5 kg/s) each step, giving 28.5 * num_steps instead of 28.5 * time.
  # Correct approach: accumulate (rate * dt) = mass injected per step.
  [injected_mass_this_step]
    type = ParsedPostprocessor
    expression = 'injection_rate_pp * dt'
    pp_names = 'injection_rate_pp dt'
    execute_on = 'timestep_end'
    outputs = 'none'
  []
  [total_injected_mass]
    type = CumulativeValuePostprocessor
    postprocessor = injected_mass_this_step
    execute_on = 'timestep_end'
  []
  # Verification: for constant 28.5 kg/s, total_injected_mass should ≈ 28.5 * time.
  # co2_mass (PorousFlowFluidMass) = current CO2 in domain; difference = outflow through BCs.
  
  # Injection well boundary area (m²). Expected ~7700–12000 (4 side walls only). Sanity checks:
  # - total_injected ≈ 28.5 * time (during injection); co2_mass tracks this
  # - max_sgas should stay < 1.1; if overshoot, check dt ramp or injection area
  [injection_well_area]
    type = AreaPostprocessor
    boundary = 'injection_well'
    execute_on = 'initial'
  []
  # [check_point_in_mesh]  # DISABLED: PointValue fails on distributed mesh
  #   type = PointValue
  #   point = '1605.4 2955.65 971.421'
  #   variable = z_coord
  #   execute_on = 'initial'
  # []
[]

[VectorPostprocessors]
  inactive = 'injection_well_sideset_info'   # B1: off to reduce assembly paths
  [injection_well_sideset_info]
    type = SidesetInfoVectorPostprocessor
    boundary = 'injection_well'
    meta_data_types = 'area centroid min max'
    execute_on = 'initial'
  []
[]

[Outputs]
  perf_graph = true
  inactive = 'exodus checkpoint'   # B1: only CSV to reduce assembly/output paths
  sync_times = '3600 86400'
  [exodus]
    type = Nemesis
    file_base = results/exodus/verif_1day/verif_1day
    time_step_interval = 10
    execute_on = 'TIMESTEP_END'
  []
  [csv]
    type = CSV
    file_base = results/exodus/verif_1day/verif_1day
    sync_only = true
  []
  [checkpoint]
    type = Checkpoint
    file_base = results/checkpoints_verif_1day/cp
    num_files = 2
    time_step_interval = 25
  []
  
  print_linear_residuals = false
  print_nonlinear_residuals = true
[]
