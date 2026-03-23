% This file is just an illustration of how to use Dynare's macro variable language,
% the documentation is available at https://www.dynare.org/manual/the-model-file.html#macro-processing-language
% To run this file, use the following command: dynare illustration savemacro=willi.mod onlymacro
% The preprocessor then creates a new mod file called "willi.mod" which looks nice and clean
% Note that the mod file here is just an illustration and incomplete;
% with a full file you don't have to use savemacro and onlymacro, but can simply run it: dynare illustration

@#define COUNTRIES = ["US","DE","FR"]
@#define CLOSING_COUNTRY = "FR"
@#define INTERMEDIATE_INPUTS = false
@#define INTERNATIONAL_BONDS = ["US","DE"]
@#define CLOSING_CONDITION = "CARRYING_COST"
@#define MODELS = ["CORE","FLEX"]
@#define TANK = true
@#if TANK == true
  @#define HOUSEHOLDS = ["_savers","_nonsavers"]
@#else
  @#define HOUSEHOLDS = [""]
@#endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                         ENDOGENOUS VARIABLES                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
@#define ENDOVARS_j = [\\
    ,"gdp"             \\
    ,"y"               \\
    ,"c"               \\
    ]
@#if INTERMEDIATE_INPUTS == true
  @#define ENDOVARS_j = ENDOVARS_j | ["x"]
@#endif

@#define ENDOVARS_j_j = ["pi"]

@#define ENDOVARS_ij_INCLUDE_jj = [\\
    ,"b"                           \\
    ,"y"                           \\
    ,"p"                           \\
    ]

@#define ENDOVARS_ij_EXCLUDE_jj = [\\
    ,"rer"                         \\
    ,"dner"                        \\
    , "piEX"                       \\
    ]

var
% country-specific variables: single index
@#for y in ENDOVARS_j
  @#define str = y
  @#define str = str + "_" + COUNTRIES[1]
  @#for j in 2:length(COUNTRIES)
    @#define str = str + " " + y + "_" + COUNTRIES[j]
  @#endfor
  @{str}
@#endfor

% country-specific variables: double index
@#for y in ENDOVARS_j_j
  @#define str = y + "_" + COUNTRIES[1] + "_" + COUNTRIES[1]
  @#for j in 2:length(COUNTRIES)
    @#define str = str + " " + y + "_" + COUNTRIES[j] + "_" + COUNTRIES[j]
  @#endfor
  @{str}
@#endfor

% cross-country variables: include own index
@#for y in ENDOVARS_ij_INCLUDE_jj
  @#for i in 1:length(COUNTRIES)
    @#define str = y + "_" + COUNTRIES[i] + "_" + COUNTRIES[1]
    @#for j in 2:length(COUNTRIES)
      @#define str = str + " " + y + "_" + COUNTRIES[i] + "_" + COUNTRIES[j]
    @#endfor
  @{str}
  @#endfor
@#endfor

% cross-country variables: exclude own index
@#for y in ENDOVARS_ij_EXCLUDE_jj
  @#for i in 1:length(COUNTRIES)
    @#define str = ""
    @#for j in 1:length(COUNTRIES) when j != i
      @#if str != ""
        @#define str = str + " "
      @#endif
      @#define str = str + y + "_" + COUNTRIES[j] + "_" + COUNTRIES[i]
    @#endfor
  @{str}
  @#endfor
@#endfor

% other variables
  a_global
  walras_test
;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          EXOGENOUS VARIABLES                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
@#define EXOVARS_j = [\\
    ,"eps_a"          \\
    ,"eps_g"          \\
    ,"eps_r"          \\
    ]

varexo
% common country-specific shocks
@#for u in EXOVARS_j
  @#define str = u
  @#define str = str + "_" + COUNTRIES[1]
    @#for j in 2:length(COUNTRIES)
      @#define str = str + " " + u + "_" + COUNTRIES[j]
    @#endfor
  @{str}
@#endfor

% global shocks
  eps_a_global
;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              PARAMETERS                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters
@#define PARAMS_j = [\\
    ,"BETA"          \\
    ,"SIGMAC"        \\
    ,"SIGMAL"        \\
    ]
@#define PARAMS_j = PARAMS_j | ["RHOG","TARGET_G_to_GDP"]
@#define PARAMS_j = PARAMS_j | ["TARGET_PI"]
@#define PARAMS_j = PARAMS_j | ["RHOR","PSIRPI","PSIRGDP"]

@#if INTERMEDIATE_INPUTS == true
  @#define PARAMS_j = PARAMS_j | ["ALPHAX"]
@#endif
@#define PARAMS_ij_INCLUDE_jj = ["PHIB"]
@#define PARAMS_ij_EXCLUDE_jj = ["OMEGA"]

% world economy normalized to 1 (SIZE_@{CLOSING_COUNTRY} is a model-local variable)
@#define SIZE_STR = ""
@#for j in 1:length(COUNTRIES)
  @#if COUNTRIES[j] != CLOSING_COUNTRY
    @#if SIZE_STR != ""
      @#define SIZE_STR = SIZE_STR + " "
    @#endif
    @#define SIZE_STR = SIZE_STR + "SIZE_" + COUNTRIES[j]
  @#endif
@#endfor
  @{SIZE_STR}

% country-specific parameters
@#for p in PARAMS_j
  @#define str = p
  @#define str = str + "_" + COUNTRIES[1]
  @#for j in 2:length(COUNTRIES)
    @#define str = str + " " + p + "_" + COUNTRIES[j]
  @#endfor
  @{str}
@#endfor

% cross-country parameters: include own index
@#for p in PARAMS_ij_INCLUDE_jj
  @#for i in 1:length(COUNTRIES)
    @#define str = p + "_" + COUNTRIES[i] + "_" + COUNTRIES[1]
    @#for j in 2:length(COUNTRIES)
      @#define str = str + " " + p + "_" + COUNTRIES[i] + "_" + COUNTRIES[j]
    @#endfor
  @{str}
  @#endfor
@#endfor

% cross_country parameters: exclude own index
@#for p in PARAMS_ij_EXCLUDE_jj
  @#for i in 1:length(COUNTRIES)
    @#define str = ""
    @#for j in 1:length(COUNTRIES) when j != i
      @#if str != ""
        @#define str = str + " "
      @#endif
      @#define str = str + p + "_" + COUNTRIES[j] + "_" + COUNTRIES[i]
    @#endfor
  @{str}
  @#endfor
@#endfor

% global parameters
  TARGET_A_GLOBAL
  RHOA_GLOBAL
;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                            MODEL EQUATIONS                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% auxiliary parameters and variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
@#define str = "1"
@#for j in COUNTRIES when j != CLOSING_COUNTRY
  @#define str = str + " - SIZE_" + ((string) j)
@#endfor
// size of @{CLOSING_COUNTRY} country (world economy is normalized to 1)
#SIZE_@{CLOSING_COUNTRY} = @{str};

// preference weights in good aggregation
@#for j in COUNTRIES
  @#for i in COUNTRIES when i != j
#gamma_@{i}_@{j} = OMEGA_@{i}_@{j} * SIZE_@{i};
  @#endfor
  @#define str = "1"
  @#for i in COUNTRIES when i != j
    @#define str = str + " - gamma_" + ((string) i) + "_" + ((string) j)
  @#endfor
#gamma_@{j}_@{j} = @{str};
@#endfor

// labor market clearing
@#for j in COUNTRIES
#nd_@{j} = l_@{j};
@#endfor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% real exchange rate: identitites %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
@#for j in 1:length(COUNTRIES)
  @#for i in j:length(COUNTRIES)
    @#if j != i
[name='real exchange rate reciprocal identity of @{COUNTRIES[i]} currency in terms of @{COUNTRIES[j]} currency']
rer_@{COUNTRIES[i]}_@{COUNTRIES[j]} = 1/rer_@{COUNTRIES[j]}_@{COUNTRIES[i]};
    @#endif
  @#endfor
@#endfor

@#if length(COUNTRIES) > 2
  @#for j in 1:length(COUNTRIES) when COUNTRIES[j] != CLOSING_COUNTRY
    @#for i in 1:length(COUNTRIES) when (j>i && COUNTRIES[i]!=CLOSING_COUNTRY)
[name='no arbitrage in spot currency markets of @{CLOSING_COUNTRY} to @{COUNTRIES[i]} via @{COUNTRIES[j]}']
rer_@{COUNTRIES[i]}_@{CLOSING_COUNTRY} = rer_@{COUNTRIES[j]}_@{CLOSING_COUNTRY} / rer_@{COUNTRIES[j]}_@{COUNTRIES[i]};
    @#endfor
  @#endfor
@#endif

@#for j in COUNTRIES
  @#for i in COUNTRIES when i != j
[name='gross change in nominal exchange rate of @{i} currency in terms of @{j} currency']
rer_@{i}_@{j} / rer_@{i}_@{j}(-1) = dner_@{i}_@{j} * pi_@{i} / pi_@{j};
  @#endfor
@#endfor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimal domestic bond decision %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
@#for MODEL in MODELS

  @#if MODEL == "CORE"
    @#define m = ""
  @#elseif MODEL == "FLEX"
    @#define m = "_flex"
  @#endif

  @#for j in COUNTRIES
    @#define h = ""
    @#define H = ""
    @#if TANK == true
      @#define h = HOUSEHOLDS[1]
      @#define H = "_SAVERS"
    @#endif
    @#define dividedByPiPlusOne = ""
    @#if MODEL == "CORE"
      @#define dividedByPiPlusOne = " / pi_" + j + "(+1)"
    @#endif
[name='optimal domestic state-contingent bond decision in country @{j}']
1 = BETA@{H}_@{j} * ( lambda@{h}@{m}_@{j}(+1) / lambda@{h}@{m}_@{j} ) * rnom@{m}_@{j}@{dividedByPiPlusOne};
  @#endfor

@#endfor // MODELS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uncovered interest rate parity %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

@#for MODEL in MODELS

  @#if MODEL == "CORE"
    @#define m = ""
  @#elseif MODEL == "FLEX"
    @#define m = "_flex"
  @#endif

  @#for j in COUNTRIES
    @#for i in COUNTRIES when j != i
      @#define dnerijPlusOne = "dner" + "_" + i + "_" + j + "(+1)"
      @#define timesdnerijPlusOne = " * " + dnerijPlusOne
      @#if MODEL == "FLEX"
        @#define dnerijPlusONe = "rer" + m + "_" + i + "_" + j + "(+1) / rer" + m + "_" + i + "_" + j
        @#define timesdnerijPlusOne = ""
      @#endif
      @#if i in INTERNATIONAL_BONDS
[name='optimal choice of @{i}-currency bonds in country @{j}']
        @#if CLOSING_CONDITION == "CARRYING_COST"
@{dnerijPlusOne} = rnom@{m}_@{j} / rnom@{m}_@{i} * ( 1 + PHIB_@{i}_@{j} * ( rer@{m}_@{i}_@{j} * b@{m}_@{i}_@{j} - steady_state(rer@{m}_@{i}_@{j}) * steady_state(b@{m}_@{i}_@{j}) ) );
        @#elseif CLOSING_CONDITION == "DEBT_ELASTIC"
1 = rnom@{m}_@{i}_@{j} / rnom@{m}_@{j}@{timesdnerijPlusOne};
        @#endif
      @#else
[name='no trade of @{i}-currency bonds in country @{j}']
b@{m}_@{i}_@{j} = 0;
      @#endif
    @#endfor
@#endfor

@#endfor // MODELS
/////////////////////
// OTHER EQUATIONS //
/////////////////////

end; // model end