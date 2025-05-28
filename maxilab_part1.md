Authors: Saskia Hekker, Susmita Das, Zhao Guo, Arthur Le Saux and Noi Shitrit for MESA School Leuven 2025

# Maxilab Part 1: Impact of convective boundary mixing on period spacing of red clump stars

## Section 1: Overview
### Science Goal

As a low-mass star continues to evolve beyond the RGB, it passes through the helium flash, followed by helium ignition in its core. This phase is called the *Red Clump*. There is still some physics about these stars that remains uncertain, in particular regarding the mixing happening above their convective core. To match obervations, several explanations have been proposed over the last decades:
1. Overmixing and penetrative convection.
2. Semiconvection (Schwarzschild & Härm 1969)
3. Maximal overshoot ([Constantino et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.452..123C/abstract))

These different mixing phenomena will have different impact on a observable quantity called the period spacing, which is sensitive to the near core region. It can also be used to estimate the mass of the helium core ([Montalbán et al. 2013](https://ui.adsabs.harvard.edu/abs/2013ApJ...766..118M/abstract)).
In this maxilab the objective is to reproduce (in a much more basic way for a purpose of saving time) the work of [Noll et al. (2024)](https://ui.adsabs.harvard.edu/abs/2024A%26A...683A.189N/abstract). Therefore, you should expect your results to only match qualitatively the restults of the paper.




### What you'll learn

The primary purpose of the first part of this maxilab is to get you more familiar with the physics of convective boundary mixing (CBM) and how it can be modelled in MESA. You will see how to use the simplest already implemented prescription as well as how to implement more complexe cases. In terms of MESA usage, you will:

1. Start a project from a test case
2. Change inlist controls
3. Add variables to output files (<span style="color:green">``history.data``</span>)
3. Customize the <span style="color:purple">``pgstar``</span> window
4. Implement a new prescription for convective boundary mixing and add new history columns using <span style="color:green">``run_star_extras.f90``</span>

### Using this Guide

Throughout this guide, monospaced text on a grey background, like this:
 ```linux
echo "Hello, World!"
```
represent commands you can type into your terminal and execute, while purple monospaced text like this: <span style="color:purple">``variable = x``</span> represents values that you might type into a file and <span style="color:green">``filename.dat``</span> are for filenames.

If you're new to Fortran, here is a short document with [some examples](https://jschwab.github.io/mesa-2016/fortran.html). Don't let yourself get hung up by the Fortran; quickly ask your classmates and the TAs for help!

Every task comes with a hint. However, if you have prior experience with MESA, do attempt to complete the task on your own. The complete solution is available <span style="color:red">here</span>.

## Section 2: Getting Started
Start by downloading the online repository ([here](https://github.com/arthurlesaux/mesasummerschool2025-day4-maxilab1/blob/eb99f57cbbc218aab5646e5c4cec46e6be7398a4/maxilab1.zip)) and uncompress it.
If you want to unzip the folder using the terminal, you can use:
```linux
unzip maxilab1.zip
```
 It contains a starting model <span style="color:green">``RC_start_noMloss.mod``</span>, some inlists to control the run (<span style="color:green">``inlist``</span>, <span style="color:green">``inlist_maxilab_pt1``</span> and <span style="color:green">``inlist_pgstar_custom``</span>) and everything you need to run a MESA model.

The starting model is a based on the MESA [1M_pre_ms_to_wd](https://docs.mesastar.org/en/latest/test_suite/1M_pre_ms_to_wd.html#m-pre-ms-to-wd) test suite. In order to save computing time, it as already been evolved from pre-main-sequence to after the helium flash. We thus start this lab just before the ignition of helium in the core of the star, i.e. the Red Clump. We have also desactivated a few of the features such as rotation and mass loss at the surface, as these can be neglected for the purpose of this lab.

Before starting the lab, check you can clean and compile the code:
```linux
cd maxilab1
./clean && ./mk
```
Make sure that you manage to start the run without any issues and break the run after you see pgplot using <span style="background-color:black"><span style="color:white">`Ctrl + C`</span></span>.

First, let's look at the <span style="color:green">``inlist_maxilab_pt1``</span>. As you can see the run starts by loading a model <span style="color:green">``RC_start_noMloss.mod``</span>, which corresponds to the start of the red clump, just after the helium flash occured and the star starts to burn helium in its core.
```fortran
&star_job

      load_saved_model = .true.
      load_model_filename = 'start_RC_noMloss.mod'

/ ! end of star_job namelist
```

The run stops when the helium content in the core is less than 0.01 in terms of mass fraction. This is done using the condition:
```fortran
&controls      

      xa_central_lower_limit_species(1) = 'he4'
      xa_central_lower_limit(1) = 0.01

/ ! end of controls namelist
```
When is criterion is met the run stops and a model is saved in <span style="color:green">``end_core_he_burn_noMloss.mod``</span> using:
```fortran
&star_job

      save_model_when_terminate = .true.
      save_model_filename = 'end_core_he_burn_noMloss.mod'
      required_termination_code_string = 'xa_central_lower_limit'

/ ! end of star_job namelist
```

## Section 3: Customizing pgstar
The next step is to customize the pgstar window to plot quantities relevant for the present study. The final pgstar window will contain 6 figures and room for a text summary at the bottom of the window. The six figures are :

1. Kippenhahn diagram
2. Hertsprung-Russel diagram
3. A profile showing power from nuclear reactions data
4. A mixing profile
5. A text panel displaying information on the current model (age, luminosity...)
6. History panel showing the period spacing and central He

First, we will start but adjusting the pgstar window so it can contain all the required plots.
<details>
  <summary> <span style="color:green">*Hints*</span> </summary>

</details>

<details>
  <summary> ***Solution*** </summary>
  This is an example but you are free to customize your window set up as you wish.

```fortran
&pgstar

  ! Set up grid layout

  file_white_on_black_flag = .false.
  !Grid1_file_flag = .true.

  Grid1_win_flag = .true.

  Grid1_file_interval = 5
  Grid1_file_width = 25

  Grid1_num_cols = 10
  Grid1_num_rows = 10

  Grid1_win_width = 14
  Grid1_win_aspect_ratio = 0.5
  Grid1_xleft = 0.00
  Grid1_xright = 1.00
  Grid1_ybot = 0.00
  Grid1_ytop = 1.00

  Grid1_num_plots = 6

  / ! end of pgstar namelist
```
</details>

Then, let's add the Kippennhahn, Herztprut-Russel, Power and Mixing diagrams
<details>
  <summary> <span style="color:green">*Hints*</span> </summary>
  These are predifined plots of pgstar. You can find information on the [pgstar documentation page](https://docs.mesastar.org/en/latest/reference/pgstar.html#pgstar).
</details>

<details>
  <summary> ***Solution*** </summary>
  This is an example but you are free to customize your window set up as you wish.

```fortran
&pgstar

  ! Add Kippenhahn plot

  Grid1_plot_name(1) = 'Kipp'

  Grid1_plot_row(1) = 1
  Grid1_plot_rowspan(1) = 5
  Grid1_plot_col(1) = 1
  Grid1_plot_colspan(1) = 4

  Grid1_plot_pad_left(1) = 0.05
  Grid1_plot_pad_right(1) = 0.01
  Grid1_plot_pad_top(1) = 0.04
  Grid1_plot_pad_bot(1) = 0.05
  Grid1_txt_scale_factor(1) = 0.5

  show_TRho_Profile_legend = .true.
  show_TRho_Profile_eos_regions = .true.
  Kipp_show_mixing = .true.
  Kipp_show_burn = .true.
  Kipp_show_luminosities = .false.

    ! Add HR diagram

  Grid1_plot_name(2) = 'HR'
  Grid1_plot_row(2) = 6
  Grid1_plot_rowspan(2) = 3
  Grid1_plot_col(2) = 1
  Grid1_plot_colspan(2) = 2

  Grid1_plot_pad_left(2) = 0.05
  Grid1_plot_pad_right(2) = 0.01
  Grid1_plot_pad_top(2) = 0.02
  Grid1_plot_pad_bot(2) = 0.07
  Grid1_txt_scale_factor(2) = 0.5

  ! Add Power profile plot

  Grid1_plot_name(3) = 'Power'

  Grid1_plot_row(3) = 6
  Grid1_plot_rowspan(3) = 3
  Grid1_plot_col(3) = 3
  Grid1_plot_colspan(3) = 2

  Grid1_plot_pad_left(3) = 0.05
  Grid1_plot_pad_right(3) = 0.01
  Grid1_plot_pad_top(3) = 0.02
  Grid1_plot_pad_bot(3) = 0.07
  Grid1_txt_scale_factor(3) = 0.5

    ! Add mixing profile

  Grid1_plot_name(5) = 'Mixing'
  Grid1_plot_row(5) = 1
  Grid1_plot_rowspan(5) = 4
  Grid1_plot_col(5) = 6
  Grid1_plot_colspan(5) = 4

  Grid1_plot_pad_left(5) = 0.05
  Grid1_plot_pad_right(5) = 0.05
  Grid1_plot_pad_top(5) = 0.04
  Grid1_plot_pad_bot(5) = 0.07
  Grid1_txt_scale_factor(5) = 0.5

  / ! end of pgstar namelist
```
</details>

The text summary panel gives information about the current state of the model. You can chose what information you want to be displayed. Let's start with: age, current time step, luminosity, effective temperature, radius, mass, h1 and he4 abundances in the core, luminosity associated with h1 and he4 burning.

<details>
  <summary> <span style="color:green">*Hints*</span> </summary>
  The variables displayed in the summary text panels are taken from the history_columns.list file.
</details>

<details>
  <summary> ***Solution*** </summary>
  This is an example but you are free to customize your panel as you wish.

```fortran
&pgstar
  ! Add text panel

  Grid1_plot_name(4) = 'Text_Summary1'
  Grid1_plot_row(4) = 9
  Grid1_plot_rowspan(4) = 2
  Grid1_plot_col(4) = 1
  Grid1_plot_colspan(4) = 10

  Grid1_plot_pad_left(4) = 0.00
  Grid1_plot_pad_right(4) = 0.00
  Grid1_plot_pad_top(4) = 0.00
  Grid1_plot_pad_bot(4) = 0.00
  Grid1_txt_scale_factor(4) = 0

  Text_Summary1_name(1,1) = 'model_number'
  Text_Summary1_name(2,1) = 'star_age'
  Text_Summary1_name(3,1) = 'log_dt'
  Text_Summary1_name(4,1) = 'luminosity'
  Text_Summary1_name(5,1) = 'Teff'
  Text_Summary1_name(6,1) = 'radius'
  Text_Summary1_name(7,1) = 'log_g'
  Text_Summary1_name(8,1) = 'star_mass'

  Text_Summary1_name(1,2) = 'log_cntr_T'
  Text_Summary1_name(2,2) = 'log_cntr_Rho'
  Text_Summary1_name(3,2) = 'log_center_P'
  Text_Summary1_name(4,2) = 'center h1'
  Text_Summary1_name(5,2) = 'center he4'
  Text_Summary1_name(6,2) = ''
  Text_Summary1_name(7,2) = ''
  Text_Summary1_name(8,2) = ''

  Text_Summary1_name(1,3) = 'log_Lnuc'
  Text_Summary1_name(2,3) = 'log_Lneu'
  Text_Summary1_name(3,3) = 'log_LH'
  Text_Summary1_name(4,3) = 'log_LHe'
  Text_Summary1_name(5,3) = 'num_zones'
  Text_Summary1_name(6,3) = ''
  Text_Summary1_name(7,3) = ''
  Text_Summary1_name(8,3) = ''

  Text_Summary1_name(1,4) = 'nu_max'
  Text_Summary1_name(2,4) = 'delta_nu'
  Text_Summary1_name(3,4) = 'delta_Pg'
  Text_Summary1_name(4,4) = ''
  Text_Summary1_name(5,4) = ''
  Text_Summary1_name(6,4) = ''
  Text_Summary1_name(7,4) = ''
  Text_Summary1_name(8,4) = ''

  / ! end of pgstar namelist
```
</details>

Last, we want to add a window displaying the period spacing and the core helium content as a function of time. This is done using a History panel.
<details>
  <summary> <span style="color:green">*Hints*</span> </summary>
  Information on how to set up a History panel can be found on [this page](https://docs.mesastar.org/en/latest/reference/pgstar.html#history-panels).
</details>

<details>
  <summary> ***Solution*** </summary>

```fortran
  &pgstar
  ! Add text panel

  ! Add mode history panel

  Grid1_plot_name(6) = 'History_Panels1' !
  Grid1_plot_row(6) = 5
  Grid1_plot_rowspan(6) = 4
  Grid1_plot_col(6) = 6
  Grid1_plot_colspan(6) = 5

  History_Panels1_win_flag = .true.
  History_Panels1_num_panels = 2
  History_Panels1_xaxis_name = 'model_number'
  History_Panels1_yaxis_name(1) ='delta_Pg'
  History_Panels1_other_yaxis_name(1) = ''
  History_Panels1_yaxis_name(2) ='center_he4'
  History_Panels1_other_yaxis_name(2) = ''
  History_Panels1_automatic_star_age_units = .true.

  Grid1_plot_pad_left(6) = 0.05
  Grid1_plot_pad_right(6) = 0.05
  Grid1_plot_pad_top(6) = 0.04
  Grid1_plot_pad_bot(6) = 0.07
  Grid1_txt_scale_factor(6) = 0.5
  / ! end of pgstar namelist
```
</details>

Now that your pgplot window is ready, let's look at the convective boundary mixing schemes.

## Section 4: Convective Boundary Mixing
### Overmixing
This scheme is already implemented in MESA, is it called the <span style="color:purple">``overshoot_scheme(:)``</span>. The aim is to add overmixing just above the helium burning convective core. Use the step overshoot scheme, over a distance of one pressure scale height, 1 H_p.

<details>
  <summary> <span style="color:green">*Hints*</span> </summary>
  Information on the overmixing scheme can be found on [this page](https://docs.mesastar.org/en/latest/reference/controls.html#overshooting). Becareful, the overmixing distance is set with <span style="color:purple">``overshoot_f(:)``</span>, which is not the distance from the convective boundary.

</details>

<details>
  <summary> ***Solution*** </summary>

```fortran
    &controls

    overshoot_scheme(1) = 'step'
    overshoot_zone_type(1) = 'burn_He'
    overshoot_zone_loc(1) = 'core'
    overshoot_bdy_loc(1) = 'top'
    overshoot_f(1) = 1.01
    overshoot_f0(1) = 0.01

    / ! end of controls namelist
```
</details>

Then, run the model with
```linux
./rn
```
The pgplot window will appear. On the right, in the mixing panel, you can see the two convective zones (in blue), which are the convective core and the envelope. Then, once helium starts burning, overshooting is taken into account (in white). As you can see, the size of these regions change with time. You will probably see additional convective zones appearing times to times, those are numerical artifacts that should not exist in a real star. However, this is a well know issue and you can just ignore them. It will not impact the results of the lab.

Just below the mixing panel, in the history panel is displayed the evolution of the period spacing (upper) and the abundance of helium at the center of the star (lower). The period spacing starts by increasing, reaches a maximum and then decrease until the end of the run. We will compare this quantity for all the runs at the end of the lab.
As expected, the helium center abundance decrease with time, as the star bruns helium.

Information about nuclear burning is provided on the Power panel. You can see that different nuclear reactions occurs at different locations. This will be studied in Maxilab2.

The run will stop once the helium core content is below 0.01, i.e. <span style="color:purple">``center_he4 < 0.01``</span>.

Before getting to the next steps, change the name of the <span style="color:green">``history.data``</span> file to, which is located in the <span style="color:green">``LOGS``</span> directory, to <span style="color:green">``history_1M_overmixing.data``</span>.

Now, we will run the three other models of CBM and compare all of them at the end of the lab.

### Penetrative convection
Penetrative convection is the same as overmixing except that we impose the temperature gradient to be adiabatic in the overshooting region. This is not implemented in MESA and should be done in a new subroutine <span style="color:purple">``penetrative_convection``</span> of <span style="color:green">``run_star_extras.f90``</span> using <span style="color:purple">``s% other_adjust_mlt_gradT_fraction``</span>.
You will need two extra subroutines. The first one <span style="color:purple">``eval_conv_bdy_Hp_perso``</span> is used to evaluate the pressure scale height at the convective boundary. The second one <span style="color:purple">``eval_over_bdy_params_perso``</span> is for evaluating other parameters such as cell index k, radius r and diffusion coefficients D. These are already implemented in <span style="color:green">``run_star_extras.f90``</span>.

<details>
  <summary> <span style="color:green">*Hints*</span> </summary>
  The overshooting parameters should be kept identical as in the overmixing case. You need to add a logical parameter in the inlist to activate or deactivate the penetrative convection scheme.
</details>

<details>
  <summary> ***Solution*** </summary>

For the inlist:
```fortran
    &controls

    x_logical_ctrl(1) = .true.

    / ! end of controls namelist
```

In <span style="color:green">``run_star_extras.f90``</span>:
```fortran
    subroutine penetrative_convection(id, ierr)
          integer, intent(in) :: id
          integer, intent(out) :: ierr
          type(star_info), pointer :: s
          integer :: i, k, k_ob, first_model_number
          real(dp) :: Hp_cb, f0, r_ob, huh, f, r, dr, factor


          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

          if (.not. s% x_logical_ctrl(1)) return
          if (s% num_conv_boundaries == 0) return
          if (.not. s% top_conv_bdy(1)) return ! no core


          f = s%overshoot_f(1)
          f0 = s%overshoot_f0(1)

          call eval_conv_bdy_Hp_perso(s, 1, Hp_cb, ierr)
          call eval_over_bdy_params_perso(s, 1, f0, k_ob, r_ob, huh, huh, ierr)

          do k = k_ob, 1, -1
            r = s%r(k)

            dr = r - r_ob

            if (dr < f*Hp_cb) then
               factor = 1._dp
               !write(*,*) 'factor = ', factor
            else
               factor = -1d0
            endif
            s% adjust_mlt_gradT_fraction(k) = factor
          end do

       end subroutine penetrative_convection
```
</details>

Then, run the model with
```linux
./rn
```
The pgplot window will appear. The run will stop once the helium core content is below 0.01, i.e. <span style="color:purple">``center_he4 < 0.01``</span>.

Before getting to the next steps, change the name of the <span style="color:green">``history.data``</span> file to, which is located in the <span style="color:green">``LOGS``</span> directory, to <span style="color:green">``history_1M_penetrative_convection.data``</span>.

### Semiconvection
The semiconvection scheme is used instead of overshooting. It is already implemented in MESA. It is independant of the overhshooting scheme, therefore both shoud not be used at the same time!

<details>
  <summary> <span style="color:green">*Hints*</span> </summary>
Information on the semiconvection scheme can be found on ([this page](https://docs.mesastar.org/en/latest/reference/controls.html#do-conv-premix)).
</details>

<details>
  <summary> ***Solution*** </summary>

```fortran
    do_conv_premix = .true.
```
</details>

Then, run the model with
```linux
./rn
```
The pgplot window will appear. In this run, there is no overshooting so the only mixing type appearing on the mixing panel will be convection (the blue one).

The run will stop once the helium core content is below 0.01, i.e. <span style="color:purple">``center_he4 < 0.01``</span>.

Before getting to the next steps, change the name of the <span style="color:green">``history.data``</span> file to, which is located in the <span style="color:green">``LOGS``</span> directory, to <span style="color:green">``history_1M_semiconvection.data``</span>.

### Maximal Overshoot
The maximal overshoot prescription is based on the model of ([Constantino et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.452..123C/abstract)). It is not implemented in MESA, but we can use a proxy for it with the <span style="color:purple">``predictive mixing``</span> scheme of MESA. It is independant of the overhshooting scheme, therefore both shoud not be used at the same time!

<details>
  <summary> <span style="color:green">*Hints*</span> </summary>
 Information on the predictive mixing scheme can be found on ([this page]( https://docs.mesastar.org/en/latest/reference/controls.html#predictive-mix))
</details>

<details>
  <summary> ***Solution*** </summary>

```fortran
    predictive_mix(1) = .true.
    predictive_zone_type(1) = 'any'
    predictive_zone_loc(1) = 'core'
    predictive_bdy_loc(1) = 'top'
    predictive_superad_thresh(1) = 5d-3
    predictive_avoid_reversal(1) = 'he4' ! avoid the happening of CBP
```
</details>

Then, run the model with
```linux
./rn
```
The pgplot window will appear. As for the semiconvection case, there is no overshooting so the only mixing type appearing on the mixing panel will be convection (the blue one).

The run will stop once the helium core content is below 0.01, i.e. <span style="color:purple">``center_he4 < 0.01``</span>.

Before getting to the next steps, change the name of the <span style="color:green">``history.data``</span> file to, which is located in the <span style="color:green">``LOGS``</span> directory, to <span style="color:green">``history_1M_max_overshoot.data``</span>.

In the end, the aim is to compare the period spacing evolution for each model using different convective boundary mixing prescription. You can of course do it yourself, but the TAs will also do it using a Python script.
