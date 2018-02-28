Dear All,

This is a note about the use of the scripts and the references. 

Some .m files are not scripts but functions, which are called from scripts (the first word in their body is ‘function’). What you want to execute (just press ‘Run’ under the Editor tab) are the scripts only. If you open Matlab, and navigate to the working directory where you saved the scripts, then in the pane called ‘Current folder’ you can see which files are actually script files, as have they distinctive symbols. 

The following was the order in which i demonstrated the scripts for you during the tutorial.

1. trajectory_of_L84.m 

2. Poincare_section_of_L84.m 

You can set the final time on line 24 to something much bigger, such as 1e5, and save the resulting variable XE, that you can load in the next two scripts. The script will run much longer, of course. You can save XE manually in the ’Command Window’ pane; check out the help file for the command ‘save’. In any case, I have already saved a large enough data set, which you find in the folder as XE.mat.

3. D0_of_L84.m

4. D1_of_L84.m

5. LEs_of_L84.m

6. MLE_of_L84.m

There is a switch (variable ’sw’) between two different integrators. For some reason Matlab’s adaptive step size ode45 does not treat some trajectories accurately enough even with a relative tolerance of 1e-10. The other integrator works much slower with a fixed time step of 1e-2, but does a good job; see MLE_of_L84_01.fig.

The equations of the Lorenz 84 (L84) model you can find in the book 

Tamás_Tél,_Márton_Gruiz_Chaotic_Dynamics_An_Introduction_Based_on_Classical_Mechanics.pdf

as eq. (9.15), but also in my articles saved in the folder that use this model as a demonstration device for various things. At the appearance of the eq. of L84 in my articles you will find references to other papers that study this model. The slides of the stuff i presented at the tutorial are: 03_Fractal_geometry.pdf and 04_Lyapunov_exponents_KY_relation.pdf.

Feel free to ask me questions about any of these.

Best wishes,

Tamas