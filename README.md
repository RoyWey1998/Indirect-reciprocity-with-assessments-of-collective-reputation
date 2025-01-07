# Indirect-reciprocity-with-assessments-of-collective-reputation

This repository includes custom code for Indirect reciprocity in the public goods game with collective reputations

The code folder includes all the basic code for the simulation and calculation. The 'default' files (default.m, default_unfixedk.m, ladefault.m, and ladefault_5.m) contain default parameter values for the reputation process of strategy evolution with fixed k, unfixed k, assessment criterion (\lambda) evolution, and co-evolution of lambda and strategy, respectively.

The 'main' files (reputation_main.m, rep_unfixed_main.m, lamain.m, and lamain_5.m) is the main processing code for these four projects. 

The files repevol2cr.m, larepevol2cr.m, larepevol2cr_5.m is the iteration function to obtian the cooperation rate from reputation evolution (for strategy evolving, lambda evolving, and co-evolving). They separately involve the corresponding functions of inner-group game functions (inner_game.m, lainner_game.m, lainner_game_5.m) and evaluation function of outsiders (outer_eval.m, laouter_eval.m, laouter_eval_5.m). 

The file repevo2cr_pofvrf.m returns payoff from simulations as the verification function. 

The files cr2pof.m, lacr2pof.m, lacr2pof_5.m use the cooperation rate to calculate expected payoffs, with the basic function payoff.m. 

The files tool_code.m and latool_code.m are collections of data-processing code in strategy and lambda evolutionary process. We calculate the selection-mutation equilibrium under rare mutations with section 2 in tool_code.m and section 1 in latool_code.m. Within the process, fixprob.m serves as the function calculating the fixation probabilities.

The file repevo_rec.m is the reputation process for recovery analysis, and the file recovery_prob_time.m is the collection of numerical calculation code for recovery analysis.

The game_pop_structure folder contains the raw data of possible population (strategy) compositions and game compositions of different population size and group size. Note that, this folder does not contain the game composition probability results for co-evolution process, which is too large to upload. Please run the code block in section 3 of latool_code.m to obtain them.

The parsave folder contains the custom saving functions in parfor loop. 

Please contact the author for any question about detailed information.
