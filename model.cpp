// Cale Seymour
// University of Nevada Las Vegas
// Written mostly in early 2018
// Carbon metabolism model written in C++.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <cstdio>
#include <string>
#define ARMA_NO_DEBUG

void mat_push(arma::mat & v, double value)
{
    arma::mat av(1,1);
    av.at(0) = value;
    v.insert_cols(v.n_cols, av.row(0));
}

// Class prototype for model.
class model
{
    public:

    // define transition matrices
    arma::mat::fixed<6,6> glucose_gP_biom_matrix; 
    arma::mat::fixed<6,6> glucose_gluconate_biom_matrix;
    arma::mat::fixed<5,6> gluconate_ribulose5PCO22NADPH_biom_matrix;
    arma::mat::fixed<1,6> gluconate_ribulose5PCO22NADPH_co2_matrix;
    arma::mat::fixed<3,6> gluconate_pyruvate_biom_matrix;
    arma::mat::fixed<3,6> gluconate_glyceraldehyde_biom_matrix;
    arma::mat::fixed<6,6> glucose_fructose_biom_matrix;
    arma::mat::fixed<6,10> ribulose5P_fructose_1_biom_matrix;
    arma::mat::fixed<6,10> ribulose5P_fructose_2_biom_matrix;
    arma::mat::fixed<3,5> ribulose5P_glyceraldehyde_biom_matrix;
    arma::mat::fixed<3,6> fructose_2glyceraldehyde_1_biom_matrix;
    arma::mat::fixed<3,6> fructose_2glyceraldehyde_2_biom_matrix;
    arma::mat::fixed<3,3> glyceraldehyde_pyruvate2NADH2ATP_biom_matrix;
    arma::mat::fixed<1,3> pyruvate_acCoACO2NADH_co2_matrix;
    arma::mat::fixed<2,3> pyruvate_acCoACO2NADH_biom_matrix;
    arma::mat::fixed<4,3> pyruvateCO2ATP_OAA_biom_matrix;
    arma::mat::fixed<6,2> acCoAOAA_ICIT_1_biom_matrix;
    arma::mat::fixed<6,4> acCoAOAA_ICIT_2_biom_matrix;
    arma::mat::fixed<5,6> ICIT_KGCO2NADH_biom_matrix;
    arma::mat::fixed<1,6> ICIT_KGCO2NADH_co2_matrix;
    arma::mat::fixed<2,6> ICIT_GLYOX_biom_matrix;
    arma::mat::fixed<4,6> ICIT_SUCC_biom_matrix;
    arma::mat::fixed<4,4> GLYOXAcCoA_OAA_biom_matrix;
    arma::mat::fixed<4,5> KG_OAACO2NADHATPFADH2_1_biom_matrix;
    arma::mat::fixed<4,5> KG_OAACO2NADHATPFADH2_2_biom_matrix;
    arma::mat::fixed<1,5> KG_OAACO2NADHATPFADH2_co2_matrix;
    
    // Response variables.
    arma::mat co2_ratios, glucose_produced_co2;
    
    // Define carbon flux independent/dependent variables.
    double r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15;
    double br1, br2, br3, br4, br5, br6, br7, br8;
    
    bool has_been_run; // We need to run the model.
    
    // Define precursor demand and isotopomer matrices and input values
    arma::mat precursor_demand, glucose_isotopomers, input_values;
    
    // Define some parameters.
    double model_tolerance, solution;
    unsigned int MAX_ITERS;
    
    // Constructor
    model();
    
    // Function prototypes.
    void modelinequality(arma::vec);
    void modelfun(arma::vec);
    void print(bool);
};

model::model()
/* Model class constructor.
   desc: builds input matrices and sets default values for precursor demand and
         glucose isotopomers. Don't set input values because we want people to
         reset those themselves!
*/
{
    has_been_run = false;
    MAX_ITERS = 45;
    r1 = 100;
    model_tolerance = 0.00000001;
    solution = 0;

    // Read isotopomers and precursor demand.
    precursor_demand =
    { {0.324536853,0.324536853,0.47,0.181639666},
      {9.285742524,9.285742524,5.01,1.301975945},
      {11.63053377,11.63053377,5.53,1.44933317 },
      {13.19748744,13.19748744,6.28,1.718110784},
      {5.014514993,5.014514993,3.27,1.040562101},
      {7.444496984,7.444496984,4.28,1.084785632},
      {4.894136918,4.894136918,2.54,0.602929144} };                 
    glucose_isotopomers =
    { {1,1,0,0,0,0,0},
      {1,0,1,0,0,0,0},
      {1,0,0,1,0,0,0},
      {1,0,0,0,1,0,0},
      {1,0,0,0,0,1,0},
      {1,0,0,0,0,0,1} };
    
    input_values =
    { {0.16,0.16,0.16,0.16,0.16,0.17} };
      
    glucose_gP_biom_matrix =
    { {1,0,0,0,0,0},
      {0,1,0,0,0,0},
      {0,0,1,0,0,0},
      {0,0,0,1,0,0},
      {0,0,0,0,1,0},
      {0,0,0,0,0,1} };
    glucose_gluconate_biom_matrix =
    { {1,0,0,0,0,0},
      {0,1,0,0,0,0},
      {0,0,1,0,0,0},
      {0,0,0,1,0,0},
      {0,0,0,0,1,0},
      {0,0,0,0,0,1} };
    gluconate_ribulose5PCO22NADPH_biom_matrix =
    { {0,1,0,0,0,0},
      {0,0,1,0,0,0},
      {0,0,0,1,0,0},
      {0,0,0,0,1,0},
      {0,0,0,0,0,1} };
    gluconate_ribulose5PCO22NADPH_co2_matrix =
    { {1,0,0,0,0,0} };
    gluconate_pyruvate_biom_matrix =
    { {1,0,0,0,0,0},
      {0,1,0,0,0,0},
      {0,0,1,0,0,0} };
    gluconate_glyceraldehyde_biom_matrix =
    { {0,0,0,1,0,0},
      {0,0,0,0,1,0},
      {0,0,0,0,0,1} };
    glucose_fructose_biom_matrix =
    { {1,0,0,0,0,0},
      {0,1,0,0,0,0},
      {0,0,1,0,0,0},
      {0,0,0,1,0,0},
      {0,0,0,0,1,0},
      {0,0,0,0,0,1} };
    ribulose5P_fructose_1_biom_matrix =
    { {1,0,0,0,0,0,0,0,0,0},
      {0,1,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,1,0,0,0,0},
      {0,0,1,0,0,0,0,0,0,0},
      {0,0,0,1,0,0,0,0,0,0},
      {0,0,0,0,1,0,0,0,0,0} };
    ribulose5P_fructose_2_biom_matrix =
    { {1,0,0,0,0,0,0,0,0,0},
      {0,1,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,1,0,0,0},
      {0,0,1,0,0,0,0,0,0,0},
      {0,0,0,1,0,0,0,0,0,0},
      {0,0,0,0,1,0,0,0,0,0} };
    ribulose5P_glyceraldehyde_biom_matrix =
    { {0,0,1,0,0},
      {0,0,0,1,0},
      {0,0,0,0,1} };
    fructose_2glyceraldehyde_1_biom_matrix =
    { {0,0,1,0,0,0},
      {0,1,0,0,0,0},
      {1,0,0,0,0,0} };
    fructose_2glyceraldehyde_2_biom_matrix =
    { {0,0,0,1,0,0},
      {0,0,0,0,1,0},
      {0,0,0,0,0,1} };
    glyceraldehyde_pyruvate2NADH2ATP_biom_matrix =
    { {1,0,0},
      {0,1,0},
      {0,0,1} };
    pyruvate_acCoACO2NADH_co2_matrix =
    { {1,0,0} };
    pyruvate_acCoACO2NADH_biom_matrix =
    { {0,1,0},
      {0,0,1} };
    pyruvateCO2ATP_OAA_biom_matrix =
    { {1,0,0},
      {0,1,0},
      {0,0,1},
      {0,0,0} };
    acCoAOAA_ICIT_1_biom_matrix =
    { {0,0},
      {0,0},
      {0,0},
      {1,0},
      {0,1},
      {0,0} };
    acCoAOAA_ICIT_2_biom_matrix =
    { {1,0,0,0},
      {0,1,0,0},
      {0,0,1,0},
      {0,0,0,0},
      {0,0,0,0},
      {0,0,0,1} };
    ICIT_KGCO2NADH_biom_matrix =
    { {1,0,0,0,0,0},
      {0,1,0,0,0,0},
      {0,0,1,0,0,0},
      {0,0,0,1,0,0},
      {0,0,0,0,1,0} };
    ICIT_KGCO2NADH_co2_matrix =
    { {0,0,0,0,0,1} };
    ICIT_GLYOX_biom_matrix =
    { {1,0,0,0,0,0},
      {0,1,0,0,0,0} };
    ICIT_SUCC_biom_matrix =
    { {0.0,0.0,0.0,0.0,0.5,0.5},
      {0.0,0.0,0.5,0.5,0.0,0.0},
      {0.0,0.0,0.5,0.5,0.0,0.0},
      {0.0,0.0,0.0,0.0,0.5,0.5} };
    GLYOXAcCoA_OAA_biom_matrix =
    { {1,0,0,0},
      {0,1,0,0},
      {0,0,0,1},
      {0,0,1,0} };
    KG_OAACO2NADHATPFADH2_1_biom_matrix =
    { {1,0,0,0,0},
      {0,1,0,0,0},
      {0,0,1,0,0},
      {0,0,0,1,0} };
    KG_OAACO2NADHATPFADH2_2_biom_matrix =
    { {0,0,0,1,0},
      {0,0,1,0,0},
      {0,1,0,0,0},
      {1,0,0,0,0} };
    KG_OAACO2NADHATPFADH2_co2_matrix =
    { {0,0,0,0,1} };
}

// Create a global model to use.
model mod;

void model::modelinequality(arma::vec x)
{
    arma::vec out;
    // r1 is defined elsewhere.
    
    r12 = x[0]; // Entner-Doudoroff pathway
    br1 = x[1]; // Glucose out
    r13 = x[2]; // Pentose phosphate pathway
    r14 = x[3]; // Glyoxylate shunt
    
    br2 = br1 * precursor_demand(0,0);
    br3 = br1 * precursor_demand(1,0);
    br4 = br1 * precursor_demand(2,0);
    br5 = br1 * precursor_demand(3,0);
    br6 = br1 * precursor_demand(4,0);
    br7 = br1 * precursor_demand(5,0);
    br8 = br1 * precursor_demand(6,0);
    
    r9 = r12 + r13;                     // Everything except glycolysis
    r8 = br6 + br7 - r14;               // Pyruvate > oxaloacetate
    r10 = (2.0/3.0) * (r13 - br8);      // Ribulose > Fructose
    r11 = (1.0/3.0) * (r13 - br8);      // Ribulose > glyceraldehyde
    r2 = r1 - r9 - br1;                 // Glucose > Fructose
    r3 = 2.0 * ( r2 + r10 - br2 );      // Fructose > 2 glyceraldehyde-P
    r4 = r3 + r12 + r11 - br3;          // Glyceraldehyde > pyruvate
    r5 = r4 + r12 - r8 - br4;           // Pyruvate > Acetyl CoA
    r15 = r5 - br5;                     // Acetyl CoA onward.
    r6 = r15 - r14;                     // Acetyl CoA > ICIT > alpha-Ketoglutarate
    r7 = r6 - br6;
    return;
}

void model::modelfun(arma::vec x)
{
    has_been_run = true;
    arma::vec OAA_OAA;
    arma::mat
    model_divergence;
    arma::mat::fixed<6,1>
    current_glucose_isotopomer;
    arma::mat::fixed<6,1>
    glucose_biom, glucose_gluconate,glucose_gluconate_biomass;
    arma::mat::fixed<3,1>
    gluconate_pyruvate, gluconate_pyruvate_biom,
    gluconate_glyceraldehyde, gluconate_glyceraldehyde_biom;
    arma::mat::fixed<6,1>
    glucose_fructose_biom;
    arma::mat::fixed<5,1>
    gluconate_ribulose5PCO22NADPH_biom, gluconate_ribulose5P;
    arma::mat::fixed<1,1>
    gluconate_ribulose5PCO22NADPH_co2;
    arma::mat::fixed<6,1>
    ribulose5P_fructose_1, ribulose5P_fructose_2, ribulose5P_fructose, ribulose5P_fructose_biom;
    arma::mat::fixed<3,1>
    ribulose5P_glyceraldehyde, ribulose5P_glyceraldehyde_biom;
    arma::mat::fixed<6,1>
    fructose_mixed;
    arma::mat::fixed<3,1>
    fructose_2glyceraldehyde_biom, glyceraldehyde_mixed, glyceraldehyde_pyruvate2NADH2ATP_biom, glyceraldehyde_pyruvate_1, pyruvate_mixed;
    arma::mat::fixed<2,1>
    pyruvate_acCoACO2NADH_biom, pyruvate_acCoA;
    arma::mat::fixed<1,1>
    pyruvate_acCoACO2NADH_co2;
    arma::mat::fixed<4,1>
    pyruvateCO2ATP_OAA, pyruvateCO2ATP_OAA_biom;
    arma::mat::fixed<6,1>
    acCoAOAA_ICIT;
    arma::mat::fixed<5,1>
    ICIT_KGCO2NADH_biom;
    arma::mat::fixed<1,1>
    ICIT_KGCO2NADH_co2;
    arma::mat::fixed<5,1>
    ICIT_KG;
    arma::mat::fixed<4,1>
    ICIT_OAA, ICIT_OAA_biom;
    arma::mat::fixed<2,1>
    ICIT_GLYOX;
    arma::mat::fixed<4,1>
    GLYOXAcCoA_OAA, GLYOXAcCoA_OAA_biom;
    arma::mat::fixed<4,1>
    KG_OAACO2NADHATPFADH2_1_biom, KG_OAACO2NADHATPFADH2_2_biom;
    arma::mat::fixed<1,1>
    KG_OAACO2NADHATPFADH2_co2;
    
    unsigned int niters;
    
    mod.modelinequality(x);    
    glucose_produced_co2.clear();
    for (unsigned int i = 0; i < glucose_isotopomers.n_cols; i++)
    {
        current_glucose_isotopomer = glucose_isotopomers.col(i);
        // glucose > glucose-P
        glucose_biom = (glucose_gP_biom_matrix * current_glucose_isotopomer) * r1;
        
        // glucose > gluconate
        glucose_gluconate = (glucose_gluconate_biom_matrix * current_glucose_isotopomer);
        glucose_gluconate_biomass = glucose_gluconate * r9;
        
        // gluconate > pyruvate + glyceraldehyde (Entner-Doudoroff)
        gluconate_pyruvate = (gluconate_pyruvate_biom_matrix * glucose_gluconate);
        gluconate_pyruvate_biom = gluconate_pyruvate * r12;
        
        gluconate_glyceraldehyde = (gluconate_glyceraldehyde_biom_matrix * glucose_gluconate);
        gluconate_glyceraldehyde_biom = gluconate_glyceraldehyde * r12;
        
        // glucose > fructose (Glycolysis)
        glucose_fructose_biom = (glucose_fructose_biom_matrix * current_glucose_isotopomer) * r2;
        
        // gluconate > ribulose5P + CO2 + 2 NADPH (Pentose Phosphate)
        gluconate_ribulose5PCO22NADPH_biom = (gluconate_ribulose5PCO22NADPH_biom_matrix * glucose_gluconate) * r13;
        gluconate_ribulose5P = (gluconate_ribulose5PCO22NADPH_biom_matrix * current_glucose_isotopomer);
        gluconate_ribulose5PCO22NADPH_co2 = (gluconate_ribulose5PCO22NADPH_co2_matrix * glucose_gluconate) * r13;

        // 3 ribulose > 2 fructose + glyceraldehyde
        ribulose5P_fructose_1 = (ribulose5P_fructose_1_biom_matrix * arma::join_cols(gluconate_ribulose5P, gluconate_ribulose5P));
        ribulose5P_fructose_2 = (ribulose5P_fructose_2_biom_matrix * arma::join_cols(gluconate_ribulose5P, gluconate_ribulose5P));
        
        // Get the average of the two fructose products, use that
        ribulose5P_fructose = (ribulose5P_fructose_1 + ribulose5P_fructose_2) / 2.0;
        ribulose5P_fructose_biom = ribulose5P_fructose * r10;
        
        ribulose5P_glyceraldehyde = (ribulose5P_glyceraldehyde_biom_matrix * gluconate_ribulose5P);
        ribulose5P_glyceraldehyde_biom = ribulose5P_glyceraldehyde * r11;
        
        // Fructose mixing pool
        if ((r2 + r10) != 0)
        {
            fructose_mixed = (glucose_fructose_biom + ribulose5P_fructose_biom) / (r2 + r10);
        }
        else
        {
            fructose_mixed = {{0,0,0,0,0,0}};
        }
        
        // fructose > 2glyceraldehyde
        fructose_2glyceraldehyde_biom = (((fructose_2glyceraldehyde_1_biom_matrix * fructose_mixed) + (fructose_2glyceraldehyde_2_biom_matrix * fructose_mixed)) * (1.0/2.0)) * r3;
        
        // glyceraldehyde mixing pool
        if ((r3 + r11 + r12) != 0)
        {
            glyceraldehyde_mixed = (fructose_2glyceraldehyde_biom + ribulose5P_glyceraldehyde_biom + gluconate_glyceraldehyde_biom) / (r3 + r11 + r12);
        }
        else
        {
            glyceraldehyde_mixed = {{0.0,0.0,0.0}};
        }
                    
        // glyceraldehyde > pyruvate + 2 NADH + 2 ATP
            glyceraldehyde_pyruvate2NADH2ATP_biom = (glyceraldehyde_pyruvate2NADH2ATP_biom_matrix * glyceraldehyde_mixed) * r4;
            glyceraldehyde_pyruvate_1 = (glyceraldehyde_pyruvate2NADH2ATP_biom_matrix * glyceraldehyde_mixed);
        
        // pyruvate mixing pool (Target molecule for TCA cycle)
        if ((r4 + r12) != 0)
        {
            pyruvate_mixed = (glyceraldehyde_pyruvate2NADH2ATP_biom + gluconate_pyruvate_biom) / (r4 + r12);
        } else {
            pyruvate_mixed = {{0.0,0.0,0.0}};
        }
        // pyruvate > acCoA + CO2 + NADH
            pyruvate_acCoACO2NADH_biom = (pyruvate_acCoACO2NADH_biom_matrix * pyruvate_mixed) * r5;
            pyruvate_acCoA = (pyruvate_acCoACO2NADH_biom_matrix * pyruvate_mixed);
            pyruvate_acCoACO2NADH_co2 = (pyruvate_acCoACO2NADH_co2_matrix * pyruvate_mixed) * r5;
        
        // pyruvate + CO2 + ATP > OAA
            pyruvateCO2ATP_OAA = (pyruvateCO2ATP_OAA_biom_matrix * pyruvate_mixed);
            pyruvateCO2ATP_OAA[3] = 0.0; //c(0,0,0,1) * 0
            pyruvateCO2ATP_OAA_biom = pyruvateCO2ATP_OAA * r8;
            
        // Run the Krebs Cycle and glyoxylate shunt
        niters = 0;
        OAA_OAA = {2,2,2,2};
        do
        {
            model_divergence = OAA_OAA;
            // acCoA + OAA > ICIT
            acCoAOAA_ICIT = (acCoAOAA_ICIT_1_biom_matrix * pyruvate_acCoA) + (acCoAOAA_ICIT_2_biom_matrix * OAA_OAA);
            
            // ICIT > KG + CO2 + NADH
            ICIT_KGCO2NADH_biom = (ICIT_KGCO2NADH_biom_matrix * acCoAOAA_ICIT) * r6;
            ICIT_KGCO2NADH_co2 = (ICIT_KGCO2NADH_co2_matrix * acCoAOAA_ICIT) * r6;
            ICIT_KG = (ICIT_KGCO2NADH_biom_matrix * acCoAOAA_ICIT);
            
            // Glyoxylate shunt
            ICIT_OAA = (ICIT_SUCC_biom_matrix * acCoAOAA_ICIT); // As far as I can tell, there is no difference in carbon arrangements between Oxaloacetate and succinate. I'll need to check on this.
            ICIT_OAA_biom = (ICIT_OAA * r14) / 2.0;
            ICIT_GLYOX = (ICIT_GLYOX_biom_matrix * acCoAOAA_ICIT);
            
            GLYOXAcCoA_OAA = (GLYOXAcCoA_OAA_biom_matrix * join_cols(ICIT_GLYOX,pyruvate_acCoA));
            GLYOXAcCoA_OAA_biom = (GLYOXAcCoA_OAA * r14) / 2.0;
            
            // KG > OAA + CO2 + NADH + ATP + FADH2
            KG_OAACO2NADHATPFADH2_1_biom = (KG_OAACO2NADHATPFADH2_1_biom_matrix * (ICIT_KG)/2.0) * r7;
            KG_OAACO2NADHATPFADH2_2_biom = (KG_OAACO2NADHATPFADH2_2_biom_matrix * (ICIT_KG)/2.0) * r7;
            KG_OAACO2NADHATPFADH2_co2 = (KG_OAACO2NADHATPFADH2_co2_matrix * ICIT_KG) * r7;
            
            // OAA
            if (r8 + r7 + r14 != 0)
            {
                OAA_OAA = (pyruvateCO2ATP_OAA_biom + KG_OAACO2NADHATPFADH2_1_biom + KG_OAACO2NADHATPFADH2_2_biom + ICIT_OAA_biom + GLYOXAcCoA_OAA_biom) / (r8 + r7 + r14);
            }
            else
            {
                OAA_OAA = (pyruvateCO2ATP_OAA_biom + KG_OAACO2NADHATPFADH2_1_biom + KG_OAACO2NADHATPFADH2_2_biom + ICIT_OAA_biom + GLYOXAcCoA_OAA_biom) * 0;
            }
            model_divergence = (OAA_OAA % OAA_OAA) - (model_divergence % model_divergence);
            niters++;
        }
        while (any(vectorise(model_divergence % model_divergence) > model_tolerance)
                 and niters < MAX_ITERS);
            // Add up the CO2 from the loop
            mat_push(glucose_produced_co2, as_scalar(gluconate_ribulose5PCO22NADPH_co2) + as_scalar(pyruvate_acCoACO2NADH_co2) + as_scalar(ICIT_KGCO2NADH_co2) + as_scalar(KG_OAACO2NADHATPFADH2_co2));
            //Rcpp::Rcout << niters << "|";
    }
    // Get isotopomer ratios
    if (glucose_produced_co2[0] != 0)
    {
        co2_ratios =
        { {glucose_produced_co2(1)/glucose_produced_co2(0),
           glucose_produced_co2(2)/glucose_produced_co2(0),
           glucose_produced_co2(3)/glucose_produced_co2(0),
           glucose_produced_co2(4)/glucose_produced_co2(0),
           glucose_produced_co2(5)/glucose_produced_co2(0),
           glucose_produced_co2(6)/glucose_produced_co2(0)} };
    } else {
        co2_ratios =
        { {0,0,0,0,0,0} };
    }
    
    //log absolute percent error
    //solution = accu(log(((co2_ratios+(1/1000000)) / input_values)));
    
    //ssq
    solution = accu((co2_ratios - input_values) % (co2_ratios - input_values));
    return;
}

void model::print(bool verbose)
{
    if (!has_been_run)
    {
        Rcpp::Rcout << "Please run the model before trying to print.\n";
        Rcpp::Rcout.flush();
        return;
    }
    double energy_production = r4 * (4.5) + (r5 * 2.5) + (r6 * 2.5) + (r7 * 7.5) + (r8 * -1.0) + (r9 * 5.0) - (r1) - (r3 / 2.0);
    double carbon_use_efficiency = (600.0 - glucose_produced_co2(0)) / 600.0;
    
    if (!verbose)
    {
        Rcpp::Rcout << r12 << "\t"
        << br1 << "\t"
        << r13 << "\t"
        << r14 << "\t"
        << solution << "\t"
        << carbon_use_efficiency << "\t"
        << energy_production << "\t"
        << input_values[0] << "\t"
        << input_values[1] << "\t"
        << input_values[2] << "\t"
        << input_values[3] << "\t"
        << input_values[4] << "\t"
        << input_values[5] << "\t"
        << co2_ratios[0] << "\t"
        << co2_ratios[1] << "\t"
        << co2_ratios[2] << "\t"
        << co2_ratios[3] << "\t"
        << co2_ratios[4] << "\t"
        << co2_ratios[5] << "\t"
        << r1 << "\t"
        << r2 << "\t"
        << r3 << "\t"
        << r4 << "\t"
        << r5 << "\t"
        << r6 << "\t"
        << r7 << "\t"
        << r8 << "\t"
        << r9 << "\t"
        << r10 << "\t"
        << r11 << "\t"
        << r12 << "\t"
        << r13 << "\t"
        << r14 << "\t"
        << r15 << "\t"
        << br1 << "\t"
        << br2 << "\t"
        << br3 << "\t"
        << br4 << "\t"
        << br5 << "\t"
        << br6 << "\t"
        << br7 << "\t"
        << br8 << "\t\n";
    } else {
        Rcpp::Rcout << "r12: " << r12 << "\t"
        << "br1: " << br1 << "\t"
        << "r13: " << r13 << "\t"
        << "r14: " << r14 << "\n"
        << "error: " << solution << "\n\n"
        << "CUE: " << carbon_use_efficiency << "\n"
        << "ENERGY: " << energy_production << "\n\n"
        << "OBSERVED:\n\t"
        << "[C1/CU] " << input_values[0] << "\n\t"
        << "[C2/CU] " << input_values[1] << "\n\t"
        << "[C3/CU] " << input_values[2] << "\n\t"
        << "[C4/CU] " << input_values[3] << "\n\t"
        << "[C5/CU] " << input_values[4] << "\n\t"
        << "[C6/CU] " << input_values[5] << "\n\t"
        << "ESTIMATED:\n\t"
        << "[C1/CU] " << co2_ratios[0] << "\n\t"
        << "[C2/CU] " << co2_ratios[1] << "\n\t"
        << "[C3/CU] " << co2_ratios[2] << "\n\t"
        << "[C4/CU] " << co2_ratios[3] << "\n\t"
        << "[C5/CU] " << co2_ratios[4] << "\n\t"
        << "[C6/CU] " << co2_ratios[5] << "\n\n"
        << "RXN FLUX:"
        << "\n\t r1|glucose in: " << r1
        << "\n\t r2|glucose > fructose (CLASSICAL GLYCOLYSIS FLUX): " << r2
        << "\n\t r3|fructose > glyceraldehyde: " << r3
        << "\n\t r4|glyceraldehyde > pyruvate: " << r4
        << "\n\t r5|pyruvate > acetylCoA: " << r5
        << "\n\t r6|isocitrate > alpha-Ketoglutarate: " << r6
        << "\n\t r7|malate > Oxaloacetate: " << r7
        << "\n\t r8|pyruvate > Oxaloacetate: " << r8
        << "\n\t r9|glucose > gluconate: " << r9
        << "\n\tr10|ribulose > fructose: " << r10
        << "\n\tr11|ribulose > glyceraldehyde: " << r11
        << "\n\tr12|gluconate > KDPG (ENTNER DOUDOROFF FLUX)*: " << r12
        << "\n\tr13|gluconate > ribulose (PENTOSE PHOSPHATE FLUX)*: " << r13
        << "\n\tr14|isocitrate + acetylCoA > 2 OAA (GLYOXYLATE FLUX)*: " << r14
        << "\n\tr15|acetylCoA > TCA and GLYOXYLATE PATHWAY: " << r15
        << "\n\tbr1|glucose out*: " << br1
        << "\n\tbr2|fructose out: " << br2
        << "\n\tbr3|glyceraldehyde out: " << br3
        << "\n\tbr4|pyruvate out: " << br4
        << "\n\tbr5|acetylCoA out: " << br5
        << "\n\tbr6|alpha-Ketoglutarate out: " << br6
        << "\n\tbr7|Oxaloacetate out: " << br7
        << "\n\tbr8|ribulose out: " << br8 << "\n";
    }
    Rcpp::Rcout.flush();
    return;
}

// [[Rcpp::export]]
double modelFunctionCpp(arma::vec x)
{
/* modelFunctionCpp:
   desc: Model function. Runs the model on the globalm model object. Accessor
         for R. (use through R, or from the global scope.)
   param: A numeric vector, representing the flux of multiple reactions.
   out: A double, representing the heuristic
*/
    mod.modelfun(x);
    //Rcpp::Rcout << ".";
    return mod.solution;
}

double invModelFunctionCpp(arma::vec x)
{
/* modelFunctionCpp:
   desc: Model function. Runs the model on the globalm model object. Accessor
         for R. (use through R, or from the global scope.) (This one is for
         maximization functions)
   param: A numeric vector, representing the flux of multiple reactions.
   out: A double, representing the heuristic
*/
    mod.modelfun(x);
    return 1/mod.solution;
}

// [[Rcpp::export]]
arma::vec modelInequalityCpp(arma::vec x)
{
/* modelInequalityCpp:
   desc: Model inequality. Calculates the flux of the model based on the given
         parameters and returns a vector of fluxes of irreversible reactions.
         (perhaps reconsider this and evaluate it as a boolean?)
   param: A numeric vector, representing the flux of multiple reactions.
   out: A numeric vector, representing the value of each of the > 0 reactions.
*/  mod.modelinequality(x);
    arma::vec out = {mod.r2, mod.r4, mod.r6, mod.r8, mod.r15, mod.br1, mod.r1 - mod.r9, mod.r1 - mod.r2};
    return out;
}

// [[Rcpp::export]]
void modelPrint(bool verbose = false)
{
    mod.print(verbose);
    return;
}

// [[Rcpp::export]]
void setModelParam(std::string param, arma::mat val)
{
/* setParam:
   desc: Set a parameter within the model. Must be one of the setable params.
         Accessor from R.
   param: A string, representing the name of the param and a numeric matrix,
          representing the value to set. (Can pass a scalar as a 1x1 matrix)
   out: none.
*/
    if (param == "glucose_isotopomers")
    {
        mod.glucose_isotopomers = val;
        return;
    }
    
    if (param == "precursor_demand")
    {
        mod.precursor_demand = val;
        return;
    }
    
    if (param == "input_values")
    {
        mod.input_values = val;
        return;
    }
    
    if (param == "model_tolerance")
    {
        mod.model_tolerance = as_scalar(val);
        return;
    }
    
    if (param == "model_tolerance")
    {
        mod.MAX_ITERS = static_cast<unsigned int>(as_scalar(val));
        return;
    }
    
    Rcpp::Rcout << "Argument1 (parameter) not found. Please select from:\n"
                << "glucose_isotopomers\n"
                << "precursor_demand\n"
                << "input_values\n"
                << "model_iters_max (default 45)\n"
                << "glucose in (default 100)\n"
                << "model_tolerance (default 0.00001)\n";
    Rcpp::Rcout.flush();
    return;
}