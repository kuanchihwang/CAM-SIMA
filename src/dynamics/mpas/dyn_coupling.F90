module dyn_coupling
    ! Modules from CAM-SIMA.
    use cam_abortutils, only: check_allocate, endrun
    use cam_constituents, only: const_is_water_species, num_advected
    use cam_thermo, only: cam_thermo_update
    use dyn_comp, only: dyn_debug_print, reverse, mpas_dynamical_core, &
        ncells_solve
    use dynconst, only: constant_cpd => cpair, constant_g => gravit, constant_p0 => pref, &
                        constant_rd => rair, constant_rv => rh2o
    use runtime_obj, only: cam_runtime_opts
    use vert_coord, only: pver, pverp

    ! Modules from CCPP.
    use cam_ccpp_cap, only: cam_constituents_array, cam_model_const_properties
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ccpp_kinds, only: kind_phys
    use geopotential_temp, only: geopotential_temp_run
    use physics_types, only: cpairv, rairv, zvirv, &
                             lagrangian_vertical, phys_state, phys_tend
    use static_energy, only: update_dry_static_energy_run

    ! Modules from CESM Share.
    use shr_kind_mod, only: kind_cx => shr_kind_cx, kind_r8 => shr_kind_r8

    implicit none

    private
    ! Provide APIs required by CAM-SIMA.
    public :: dynamics_to_physics_coupling
    public :: physics_to_dynamics_coupling
contains
    !> This subroutine implements one-way coupling from the dynamics output states to the physics input states.
    !> The other coupling direction is implemented by its counterpart, `physics_to_dynamics_coupling`.
    !> (KCW, 2024-07-31)
    subroutine dynamics_to_physics_coupling()
        character(*), parameter :: subname = 'dyn_coupling::dynamics_to_physics_coupling'
        integer :: column_index
        integer, pointer :: index_qv
        real(kind_phys), pointer :: constituents(:, :, :)
        ! Variable name suffixes have the following meanings:
        ! `*_col`: Variable is of each column.
        ! `*_int`: Variable is at layer interfaces.
        ! `*_mid`: Variable is at layer midpoints.
        real(kind_r8), allocatable :: pd_int_col(:), &       ! Dry hydrostatic air pressure (Pa).
                                      pd_mid_col(:), &       ! Dry hydrostatic air pressure (Pa).
                                      p_int_col(:), &        ! Full hydrostatic air pressure (Pa).
                                      p_mid_col(:), &        ! Full hydrostatic air pressure (Pa).
                                      z_int_col(:)           ! Geometric height (m).
        real(kind_r8), allocatable :: dpd_col(:), &          ! Dry air pressure difference (Pa) between layer interfaces.
                                      dp_col(:), &           ! Full air pressure difference (Pa) between layer interfaces.
                                      dz_col(:)              ! Geometric height difference (m) between layer interfaces.
        real(kind_r8), allocatable :: qv_mid_col(:), &       ! Water vapor mixing ratio (kg kg-1).
                                      sigma_all_q_mid_col(:) ! Summation of all water mixing ratio (kg kg-1).
        real(kind_r8), allocatable :: rhod_mid_col(:), &     ! Dry air density (kg m-3).
                                      rho_mid_col(:)         ! Full air density (kg m-3).
        real(kind_r8), allocatable :: t_mid_col(:), &        ! Temperature (K).
                                      tm_mid_col(:), &       ! Modified "moist" temperature (K).
                                      tv_mid_col(:)          ! Virtual temperature (K).
        real(kind_r8), allocatable :: u_mid_col(:), &        ! Eastward wind velocity (m s-1).
                                      v_mid_col(:), &        ! Northward wind velocity (m s-1).
                                      omega_mid_col(:)       ! Vertical wind velocity (Pa s-1).
        real(kind_r8), pointer :: exner(:, :)
        real(kind_r8), pointer :: rho_zz(:, :)
        real(kind_r8), pointer :: scalars(:, :, :)
        real(kind_r8), pointer :: theta_m(:, :)
        real(kind_r8), pointer :: ucellzonal(:, :), ucellmeridional(:, :), w(:, :)
        real(kind_r8), pointer :: zgrid(:, :)
        real(kind_r8), pointer :: zz(:, :)

        call init_shared_variable()

        ! Set variables in the `physics_state` derived type column by column.
        ! This way, peak memory usage can be reduced.
        do column_index = 1, ncells_solve
            call update_shared_variable()
            call set_physics_state_column()
        end do

        call set_physics_state_external()

        deallocate(pd_int_col, pd_mid_col, p_int_col, p_mid_col, z_int_col)
        deallocate(dpd_col, dp_col, dz_col)
        deallocate(qv_mid_col, sigma_all_q_mid_col)
        deallocate(rhod_mid_col, rho_mid_col)
        deallocate(t_mid_col, tm_mid_col, tv_mid_col)
        deallocate(u_mid_col, v_mid_col, omega_mid_col)

        nullify(index_qv)
        nullify(constituents)
        nullify(exner)
        nullify(rho_zz)
        nullify(scalars)
        nullify(theta_m)
        nullify(ucellzonal, ucellmeridional, w)
        nullify(zgrid)
        nullify(zz)
    contains
        !> Initialize variables that are shared and repeatedly used by the `update_shared_variable` and
        !> `set_physics_state_column` internal subroutines.
        !> (KCW, 2024-07-20)
        subroutine init_shared_variable()
            character(*), parameter :: subname = 'dyn_coupling::dynamics_to_physics_coupling::init_shared_variable'
            integer :: ierr

            nullify(index_qv)
            nullify(constituents)
            nullify(exner)
            nullify(rho_zz)
            nullify(scalars)
            nullify(theta_m)
            nullify(ucellzonal, ucellmeridional, w)
            nullify(zgrid)
            nullify(zz)

            allocate(pd_int_col(pverp), pd_mid_col(pver), p_int_col(pverp), p_mid_col(pver), z_int_col(pverp), stat=ierr)
            call check_allocate(ierr, subname, &
                'pd_int_col(pverp), pd_mid_col(pver), p_int_col(pverp), p_mid_col(pver), z_int_col(pverp)', &
                'dyn_coupling', __LINE__)

            allocate(dpd_col(pver), dp_col(pver), dz_col(pver), stat=ierr)
            call check_allocate(ierr, subname, &
                'dpd_col(pver), dp_col(pver), dz_col(pver)', &
                'dyn_coupling', __LINE__)

            allocate(qv_mid_col(pver), sigma_all_q_mid_col(pver), stat=ierr)
            call check_allocate(ierr, subname, &
                'qv_mid_col(pver), sigma_all_q_mid_col(pver)', &
                'dyn_coupling', __LINE__)

            allocate(rhod_mid_col(pver), rho_mid_col(pver), stat=ierr)
            call check_allocate(ierr, subname, &
                'rhod_mid_col(pver), rho_mid_col(pver)', &
                'dyn_coupling', __LINE__)

            allocate(t_mid_col(pver), tm_mid_col(pver), tv_mid_col(pver), stat=ierr)
            call check_allocate(ierr, subname, &
                't_mid_col(pver), tm_mid_col(pver), tv_mid_col(pver)', &
                'dyn_coupling', __LINE__)

            allocate(u_mid_col(pver), v_mid_col(pver), omega_mid_col(pver), stat=ierr)
            call check_allocate(ierr, subname, &
                'u_mid_col(pver), v_mid_col(pver), omega_mid_col(pver)', &
                'dyn_coupling', __LINE__)

            constituents => cam_constituents_array()

            if (.not. associated(constituents)) then
                call endrun('Failed to find variable "constituents"', subname, __LINE__)
            end if

            call mpas_dynamical_core % get_variable_pointer(index_qv, 'dim', 'index_qv')
            call mpas_dynamical_core % get_variable_pointer(exner, 'diag', 'exner')
            call mpas_dynamical_core % get_variable_pointer(rho_zz, 'state', 'rho_zz', time_level=1)
            call mpas_dynamical_core % get_variable_pointer(scalars, 'state', 'scalars', time_level=1)
            call mpas_dynamical_core % get_variable_pointer(theta_m, 'state', 'theta_m', time_level=1)
            call mpas_dynamical_core % get_variable_pointer(ucellzonal, 'diag', 'uReconstructZonal')
            call mpas_dynamical_core % get_variable_pointer(ucellmeridional, 'diag', 'uReconstructMeridional')
            call mpas_dynamical_core % get_variable_pointer(w, 'state', 'w', time_level=1)
            call mpas_dynamical_core % get_variable_pointer(zgrid, 'mesh', 'zgrid')
            call mpas_dynamical_core % get_variable_pointer(zz, 'mesh', 'zz')
        end subroutine init_shared_variable

        !> Update variables for the specific column, indicated by `column_index`. This subroutine and `set_physics_state_column`
        !> should be called in pairs.
        !> (KCW, 2024-07-30)
        subroutine update_shared_variable()
            character(*), parameter :: subname = 'dyn_coupling::dynamics_to_physics_coupling::update_shared_variable'
            integer :: i, j, k

            i = column_index

            ! The summation term of equation 5 in doi:10.1029/2017MS001257.
            sigma_all_q_mid_col(:) = 1.0_kind_r8

            do j = 1, num_advected
                if (const_is_water_species(j)) then
                    sigma_all_q_mid_col(:) = &
                        sigma_all_q_mid_col(:) + scalars(mpas_dynamical_core % map_mpas_scalar_index(j), :, i)
                end if
            end do

            ! Compute thermodynamic variables.

            ! By definition.
            z_int_col(:) = zgrid(:, i)
            dz_col(:) = z_int_col(2:pverp) - z_int_col(1:pver)
            qv_mid_col(:) = scalars(index_qv, :, i)
            rhod_mid_col(:) = rho_zz(:, i) * zz(:, i)

            ! Equation 5 in doi:10.1029/2017MS001257.
            rho_mid_col(:) = rhod_mid_col(:) * sigma_all_q_mid_col(:)

            ! Hydrostatic equation.
            dpd_col(:) = -rhod_mid_col(:) * constant_g * dz_col(:)
            dp_col(:) = -rho_mid_col(:) * constant_g * dz_col(:)

            ! By definition of Exner function. Also see below.
            tm_mid_col(:) = theta_m(:, i) * exner(:, i)

            ! Paragraph below equation 2.7 in doi:10.5065/1DFH-6P97.
            ! Paragraph below equation 2 in doi:10.1175/MWR-D-11-00215.1.
            t_mid_col(:) = tm_mid_col(:) / &
                (1.0_kind_r8 + constant_rv / constant_rd * qv_mid_col(:))

            ! Equation 16 in doi:10.1029/2017MS001257.
            ! The numerator terms are just `tm_mid_col` here (i.e., modified "moist" temperature).
            tv_mid_col(:) = tm_mid_col(:) / sigma_all_q_mid_col(:)

            ! Hydrostatic equation with equation of state plugged in and arranging for pressure.
            pd_mid_col(:) = -constant_rd * t_mid_col(:) * dpd_col(:) / (constant_g * dz_col(:))
            p_mid_col(:) = -constant_rd * tv_mid_col(:) * dp_col(:) / (constant_g * dz_col(:))

            ! Assume no water at top of model.
            pd_int_col(pverp) = p_mid_col(pver) + 0.5_kind_r8 * dp_col(pver)
            p_int_col(pverp) = pd_int_col(pverp)

            ! Integrate downward.
            do k = pver, 1, -1
                pd_int_col(k) = pd_int_col(k + 1) - dpd_col(k)
                p_int_col(k) = p_int_col(k + 1) - dp_col(k)
            end do

            ! Compute momentum variables.

            ! By definition.
            u_mid_col(:) = ucellzonal(:, i)
            v_mid_col(:) = ucellmeridional(:, i)
            omega_mid_col(:) = -rhod_mid_col(:) * constant_g * 0.5_kind_r8 * (w(1:pver, i) + w(2:pverp, i))
        end subroutine update_shared_variable

        !> Set variables for the specific column, indicated by `column_index`, in the `physics_state` derived type.
        !> This subroutine and `update_shared_variable` should be called in pairs.
        !> (KCW, 2024-07-30)
        subroutine set_physics_state_column()
            character(*), parameter :: subname = 'dyn_coupling::dynamics_to_physics_coupling::set_physics_state_column'
            integer :: i, j

            i = column_index

            ! Vertical index order is reversed between CAM-SIMA and MPAS.
            ! Always call `reverse` when assigning anything in the `physics_state` derived type.

            phys_state % u(i, :) = reverse(u_mid_col)
            phys_state % v(i, :) = reverse(v_mid_col)
            phys_state % omega(i, :) = reverse(omega_mid_col)

            phys_state % psdry(i) = pd_int_col(1)
            phys_state % pintdry(i, :) = reverse(pd_int_col)
            phys_state % pmiddry(i, :) = reverse(pd_mid_col)
            phys_state % pdeldry(i, :) = reverse(-dpd_col)
            phys_state % lnpintdry(i, :) = log(phys_state % pintdry(i, :))
            phys_state % lnpmiddry(i, :) = log(phys_state % pmiddry(i, :))
            phys_state % rpdeldry(i, :) = 1.0_kind_r8 / phys_state % pdeldry(i, :)

            phys_state % ps(i) = p_int_col(1)
            phys_state % pint(i, :) = reverse(p_int_col)
            phys_state % pmid(i, :) = reverse(p_mid_col)
            phys_state % pdel(i, :) = reverse(-dp_col)
            phys_state % lnpint(i, :) = log(phys_state % pint(i, :))
            phys_state % lnpmid(i, :) = log(phys_state % pmid(i, :))
            phys_state % rpdel(i, :) = 1.0_kind_r8 / phys_state % pdel(i, :)

            phys_state % t(i, :) = reverse(t_mid_col)

            ! In CAM-SIMA, constituents are held by CCPP rather than the `physics_state` derived type.
            ! `constituents` points to CCPP memory.
            ! `scalars` points to MPAS memory.
            do j = 1, num_advected
                constituents(i, :, j) = &
                    reverse(scalars(mpas_dynamical_core % map_mpas_scalar_index(j), :, i))
            end do

            ! This variable name is very misleading. It is actually the reciprocal of Exner function.
            ! It should have been named `rexner` or similar...
            ! Also, do not use `exner` from MPAS here because pressure in MPAS is non-hydrostatic.
            phys_state % exner(i, :) = (constant_p0 / phys_state % pmid(i, :)) ** (constant_rd / constant_cpd)

            phys_state % phis(i) = constant_g * z_int_col(1)
        end subroutine set_physics_state_column

        !> Set variables in the `physics_state` derived type by calling external subroutines.
        !> (KCW, 2024-07-30)
        subroutine set_physics_state_external()
            character(*), parameter :: subname = 'dyn_coupling::dynamics_to_physics_coupling::set_physics_state_external'
            character(kind_cx) :: cerr
            integer :: ierr
            type(ccpp_constituent_prop_ptr_t), pointer :: constituent_properties(:)

            nullify(constituent_properties)

            constituent_properties => cam_model_const_properties()

            if (.not. associated(constituent_properties)) then
                call endrun('Failed to find variable "constituent_properties"', subname, __LINE__)
            end if

            ! Call this subroutine to update `cpairv`, `rairv`, `zvirv`, etc.
            call cam_thermo_update( &
                constituents, phys_state % t, ncells_solve, cam_runtime_opts % update_thermodynamic_variables())

            ! Set `zi` (i.e., geopotential height at layer interfaces) and `zm` (i.e., geopotential height at layer midpoints).
            call geopotential_temp_run( &
                pver, lagrangian_vertical, pver, 1, pverp, 1, num_advected, &
                phys_state % lnpint, phys_state % pint, phys_state % pmid, phys_state % pdel, phys_state % rpdel, phys_state % t, &
                constituents(:, :, mpas_dynamical_core % map_constituent_index(index_qv)), constituents, &
                constituent_properties, rairv, constant_g, zvirv, phys_state % zi, phys_state % zm, ncells_solve, ierr, cerr)

            ! Set `dse` (i.e., dry static energy).
            call update_dry_static_energy_run( &
                pver, constant_g, phys_state % t, phys_state % zm, phys_state % phis, phys_state % dse, cpairv, ierr, cerr)

            nullify(constituent_properties)
        end subroutine set_physics_state_external
    end subroutine dynamics_to_physics_coupling

    !> This subroutine implements one-way coupling from the physics output states to the dynamics input states.
    !> The other coupling direction is implemented by its counterpart, `dynamics_to_physics_coupling`.
    !> TODO
    subroutine physics_to_dynamics_coupling()
        character(*), parameter :: subname = 'dyn_coupling::physics_to_dynamics_coupling'
    end subroutine physics_to_dynamics_coupling
end module dyn_coupling
