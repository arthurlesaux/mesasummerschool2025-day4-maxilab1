! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use num_lib, only: find0, two_piece_linear_coeffs

      implicit none

      ! these routines are called by the standard run_star check_model
      contains

      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

      end subroutine extras_controls


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_startup


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.


      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do

      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! to save a profile,
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step


      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve

      !!! Write penetrative_convection subroutine here !!!

      subroutine penetrative_convection
          integer, intent(out) :: ierr
          type(star_info), pointer :: s
      end subroutine penetrative_convection

       subroutine evaluate_cz_bdy_dq(s, cz_bdy_dq)
          type(star_info), pointer, intent(in) :: s
          real(dp), intent(out) :: cz_bdy_dq(s% nz)
          real(dp) :: cz_dq
          integer :: idx_bdy

          cz_bdy_dq = 0d0
          ! We only do it for the first limit...


          ! Test and works!
          idx_bdy = s%conv_bdy_loc(1)
          cz_dq = find0(0._dp, s% gradr(idx_bdy+1) - s%gradL(idx_bdy+1),&
            s%dq(idx_bdy),  s% gradr(idx_bdy) - s% gradL(idx_bdy))

          ! If strange value (negative or bigger than dq) then don't interpolate and take the closest
          ! Probably not the best...
          cz_bdy_dq(idx_bdy) = max(0d0, min(s% dq(idx_bdy) - cz_dq, s% dq(idx_bdy)))


       end subroutine

       ! No diff with the overshoot_utils one
       subroutine eval_conv_bdy_k_perso (s, i, k, ierr)

          type(star_info), pointer :: s
          integer, intent(in)      :: i
          integer, intent(out)     :: k
          integer, intent(out)     :: ierr

          ! Evaluate the index k of the cell containing the i'th convective
          ! boundary

          ierr = 0

          if (s%top_conv_bdy(i)) then
             k = s%conv_bdy_loc(i)
          else
             k = s%conv_bdy_loc(i) - 1
          endif

          if (k >= s%nz .OR. k < 1) then
             write(*,*) 'Invalid cell for convective boundary: i, k, nz=', i, k, s%nz
             ierr = -1
             return
          endif

          ! Finish

          return

        end subroutine eval_conv_bdy_k_perso

        !****

        subroutine eval_conv_bdy_r_perso (s, i, r, ierr)

          type(star_info), pointer :: s
          integer, intent(in)      :: i
          real(dp), intent(out)    :: r
          integer, intent(out)     :: ierr

          integer  :: k
          real(dp) :: w, cz_bdy_dq(s% nz)

          ! Evaluate the radius r at the i'th convective boundary

          ! Find the convective boundary cell

          ierr = 0

          call eval_conv_bdy_k_perso(s, i, k, ierr)
          if (ierr /= 0) return

          ! Interpolate r based on the fact that r^3 varies linearly with q
          ! across the (constant-density) cell

          call evaluate_cz_bdy_dq(s, cz_bdy_dq)

          w = cz_bdy_dq(k)/s%dq(k)

          if (w < 0._dp .OR. w > 1._dp) then
             write(*,*) 'Invalid weight for convective boundary: i, w=', i, w
             ierr = -1
             return
          end if

          associate (k_o => k, &
                     k_i => k+1)

            r = pow((      w)*s%r(k_i)*s%r(k_i)*s%r(k_i) + &
                       (1._dp-w)*s%r(k_o)*s%r(k_o)*s%r(k_o), 1._dp/3._dp)

          end associate

          ! Finish

          return

        end subroutine eval_conv_bdy_r_perso

        !****

        subroutine eval_conv_bdy_Hp_perso (s, i, Hp, ierr)

          type(star_info), pointer :: s
          integer, intent(in)      :: i
          real(dp), intent(out)    :: Hp
          integer, intent(out)     :: ierr

          integer  :: k
          real(dp) :: r
          real(dp) :: w
          real(dp) :: x0
          real(dp) :: x1
          real(dp) :: x2
          real(dp) :: x
          real(dp) :: a0
          real(dp) :: a1
          real(dp) :: a2
          real(dp) :: P
          real(dp) :: rho
          real(dp) :: r_top
          real(dp) :: r_bot
          real(dp) :: dr
          real(dp) :: cz_bdy_dq(s% nz)

          ! Evaluate the pressure scale height Hp at the i'th convective boundary

          ! Find the convective boundary cell

          ierr = 0
          call evaluate_cz_bdy_dq(s, cz_bdy_dq)

          call eval_conv_bdy_k_perso(s, i, k, ierr)
          if (ierr /= 0) return

          ! Evaluate the radius at the convective boundary

          call eval_conv_bdy_r_perso(s, i, r, ierr)
          if (ierr /= 0) return

          ! Interpolate the pressure and density at the boundary, using a
          ! quadratic fit across the boundary cell and its neighbors (the
          ! x's are fractional mass distances from the outer edge of cell
          ! k-1); then, evaluate the pressure scale height

          associate (k_o => k-1, &
                     k_m => k, &
                     k_i => k+1)

            x0 = s%dq(k_o)/2._dp
            x1 = s%dq(k_o) + s%dq(k_m)/2._dp
            x2 = s%dq(k_o) + s%dq(k_m) + s%dq(k_i)/2._dp

            x = s%dq(k_o) + cz_bdy_dq(k)

            call two_piece_linear_coeffs(x, x0, x1, x2, a0, a1, a2, ierr)
            if (ierr /= 0) return

            P = exp(a0*s%lnPeos(k_o) + a1*s%lnPeos(k_m) + a2*s%lnPeos(k_i))
            rho = exp(a0*s%lnd(k_o) + a1*s%lnd(k_m) + a2*s%lnd(k_i))

            ! Evaluate the pressure scale height

            Hp = P/(rho*s%cgrav(k_m)* &
                 (s%M_center + s%xmstar*s%conv_bdy_q(i))/(r*r))

          end associate

          ! (Possibly) limit the scale height using the size of the
          ! convection zone

          if (s%limit_overshoot_Hp_using_size_of_convection_zone) then

             ! Determine the radial extent of the convection zone (note that
             ! r_top/r_bot don't coincide exactly with the r calculated
             ! above)

             if (s%top_conv_bdy(i)) then

                if (i == 1) then
                   r_bot = s%R_center
                else
                   if (s%top_conv_bdy(i-1)) then
                      write(*,*) 'Double top boundary in overshoot; i=', i
                      ierr = -1
                      return
                   end if
                   r_bot = s%r(s%conv_bdy_loc(i-1))
                endif

                r_top = s%r(k)

             else

                r_bot = s%r(k+1)

                if (i == s%num_conv_boundaries) then
                   r_top = s%r(1)
                else
                   if (.NOT. s%top_conv_bdy(i+1)) then
                      write(*,*) 'Double bottom boundary in overshoot; i=', i
                      ierr = -1
                      return
                   endif
                   r_top = s%r(s%conv_bdy_loc(i+1))
                endif

             endif

             dr = r_top - r_bot

             ! Apply the limit

             if (s%overshoot_alpha > 0d0) then
                if (s%overshoot_alpha*Hp > dr) Hp = dr/s%overshoot_alpha
             else
                if (s%alpha_mlt(k)*Hp > dr) Hp = dr/s%mixing_length_alpha
             end if

          end if

          ! Finish

          return

        end subroutine eval_conv_bdy_Hp_perso

        !****

        subroutine eval_over_bdy_params_perso (s, i, f0, k, r, D, vc, ierr)

          type(star_info), pointer :: s
          integer, intent(in)      :: i
          real(dp), intent(in)     :: f0
          integer, intent(out)     :: k
          real(dp), intent(out)    :: r
          real(dp), intent(out)    :: D
          real(dp), intent(out)    :: vc
          integer, intent(out)     :: ierr

          integer  :: k_cb
          real(dp) :: r_cb
          real(dp) :: Hp_cb
          real(dp) :: w
          real(dp) :: lambda

          ! Evaluate parameters (cell index k, radius r, diffusion
          ! coefficients D and cdc) for the overshoot boundary associated
          ! with the i'th convective boundary

          ! Find the convective boundary cell

          ierr = 0

          call eval_conv_bdy_k_perso(s, i, k_cb, ierr)
          if (ierr /= 0) return

          ! Evaluate the radius at the convective boundary

          call eval_conv_bdy_r_perso(s, i, r_cb, ierr)
          if (ierr /= 0) return

          ! Evaluate the pressure scale height at the convective boundary

          call eval_conv_bdy_Hp_perso(s, i, Hp_cb, ierr)
          if (ierr /= 0) return

          ! Search for the overshoot boundary cell

          ierr = 0

          if (s%top_conv_bdy(i)) then

             ! Overshooting outward -- search inward

             r = r_cb - f0*Hp_cb

             if (r <= s%r(s%nz)) then

                r = s%r(s%nz)
                k = s%nz - 1

             else

                search_in_loop: do k = k_cb, s%nz-1
                   if (s%r(k+1) <= r) exit search_in_loop
                end do search_in_loop

             endif

          else

             ! Overshooting inward -- search outward

             r = r_cb + f0*Hp_cb

             if (r >=  s%r(1)) then

                r = s%r(1)
                k = 1

             else

                search_out_loop : do k = k_cb, 1, -1
                   if (s%r(k) > r) exit search_out_loop
                end do search_out_loop

             endif

          endif

          if (.NOT. (s%r(k+1) <= r .AND. s%r(k) >= r)) then
             write(*,*) 'r_ob not correctly bracketed: r(k+1), r, r(k)=', s%r(k+1), r, s%r(k)
             ierr = -1
             return
          end if

          ! Interpolate mixing parameters

          w = (s%r(k)*s%r(k)*s%r(k) - r*r*r)/ &
              (s%r(k)*s%r(k)*s%r(k) - s%r(k+1)*s%r(k+1)*s%r(k+1))

          lambda = (1._dp-w)*s%mlt_mixing_length(k) + w*s%mlt_mixing_length(k+1)

          if (s%conv_vel(k) /= 0._dp .AND. s%conv_vel(k+1) /= 0._dp) then

             ! Both faces of cell have non-zero mixing; interpolate vc between faces

             vc = (1._dp-w)*s%conv_vel(k) + w*s%conv_vel(k+1)

          elseif (s%conv_vel(k) /= 0._dp .AND. s%conv_vel(k+1) == 0._dp) then

             ! Outer face of cell has non-zero mixing; interpolate vc
             ! between this face and r_cb, assuming vc = 0 at the latter

              if(s%r(k) /= r_cb) then
                w = (s%r(k)*s%r(k)*s%r(k) - r*r*r)/ &
                 (s%r(k)*s%r(k)*s%r(k) - r_cb*r_cb*r_cb)
              else
                w = 0d0
              end if

             vc = (1._dp-w)*s%conv_vel(k)

          elseif (s%conv_vel(k) == 0._dp .AND. s%conv_vel(k+1) /= 0._dp) then

             ! Inner face of cell has non-zero mixing; interpolate vc
             ! between this face and r_cb, assuming vc = 0 at the latter

             if(s%r(k+1) /= r_cb) then
                w = (r_cb*r_cb*r_cb - r*r*r)/ &
                 (r_cb*r_cb*r_cb - s%r(k+1)*s%r(k+1)*s%r(k+1))
             else
                w = 0d0
             end if

             vc = w*s%conv_vel(k+1)

          else

             ! Neither face of cell has non-zero mixing; return

             vc = 0._dp

          endif

          ! Evaluate the diffusion coefficient

          D = vc*lambda/3._dp

          ! Finish

          ierr = 0

          return

        end subroutine eval_over_bdy_params_perso


      end module run_star_extras
