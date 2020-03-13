//  Author: Robert Scheller, Melissa Lucash

using Landis.Core;
using Landis.SpatialModeling;
using Landis.Utilities;
using System;
using Landis.Library.Climate;

namespace Landis.Extension.Succession.NECN
{
    public enum SoilName { Primary, Secondary, Tertiary, Other };
    /// <summary>
    /// </summary>
    public class SoilLayer 
    {

        private SoilName name;
        //private LayerType type;
        private double o_carbon;
        private double o_nitrogen;
        private double do_carbon;
        private double do_nitrogen;
        private double microbial_carbon;
        private double microbial_nitrogen;
        private double enzymatic_concentration;
        private double monthlyCinputs;
        private double monthlyNinputs;
        //private double decayValue;
        //private double fractionLignin;
        //private double netMineralization;
        //private double grossMineralization;


        //---------------------------------------------------------------------
        public SoilLayer(SoilName name)//, LayerType type)
        {
            this.name = name;
            //this.type = type;
            this.o_carbon = 0.0;
            this.o_nitrogen = 0.0;
            this.do_carbon = 0.0;
            this.do_nitrogen = 0.0;
            this.microbial_carbon = 0.0;
            this.microbial_nitrogen = 0.0;
            this.enzymatic_concentration = 0.0;
            this.monthlyCinputs = 0.0;

        }
        //---------------------------------------------------------------------
        /// <summary>
        /// Layer Name
        /// </summary>
        public SoilName Name
        {
            get
            {
                return name;
            }
            set
            {
                name = value;
            }
        }
        //---------------------------------------------------------------------
        /// <summary>
        /// Soil Organic Carbon
        /// </summary>
        public double Carbon
        {
            get
            {
                return o_carbon;
            }
            set
            {
                o_carbon = value;
            }
        }
        //---------------------------------------------------------------------
        /// <summary>
        /// Soil Organic Carbon
        /// </summary>
        public double MonthlyCarbonInputs
        {
            get
            {
                return monthlyCinputs;
            }
            set
            {
                monthlyCinputs = value;
            }
        }
        //---------------------------------------------------------------------
        /// <summary>
        /// Soil Organic Carbon
        /// </summary>
        public double MonthlyNitrogenInputs
        {
            get
            {
                return monthlyNinputs;
            }
            set
            {
                monthlyNinputs = value;
            }
        }
        //---------------------------------------------------------------------

        /// <summary>
        /// Soil Organic Nitrogen
        /// </summary>
        public double Nitrogen
        {
            get
            {
                return o_nitrogen;
            }
            set
            {
                o_nitrogen = value;
            }
        }
        //---------------------------------------------------------------------
        /// <summary>
        /// Dissolved Organic Carbon
        /// </summary>
        public double DOC
        {
            get
            {
                return do_carbon;
            }
            set
            {
                do_carbon = value;
            }
        }
        //---------------------------------------------------------------------

        /// <summary>
        /// Dissolved Organic Nitrogen
        /// </summary>
        public double DON
        {
            get
            {
                return do_nitrogen;
            }
            set
            {
                do_nitrogen = value;
            }
        }
        //---------------------------------------------------------------------
        /// <summary>
        /// Microbial Carbon
        /// </summary>
        public double MicrobialCarbon
        {
            get
            {
                return microbial_carbon;
            }
            set
            {
                microbial_carbon = value;
            }
        }
        //---------------------------------------------------------------------

        /// <summary>
        /// Microbial Nitrogen
        /// </summary>
        public double MicrobialNitrogen
        {
            get
            {
                return microbial_nitrogen;
            }
            set
            {
                microbial_nitrogen = value;
            }
        }
        //---------------------------------------------------------------------

        /// <summary>
        /// Enzymatic Concentration
        /// </summary>
        public double EnzymaticConcentration
        {
            get
            {
                return enzymatic_concentration;
            }
            set
            {
                enzymatic_concentration = value;
            }
        }

        // Moved to Input file for calibration
        private static double r_death = PlugIn.Parameters.MicrobialTurnoverRate;
        private static double frac = PlugIn.Parameters.FractionUnprotectedSOM;
        private static double r_ecloss = PlugIn.Parameters.EnzymeTurnoverRate;
        private static double c_use_efficiency = PlugIn.Parameters.CarbonUseEfficiency;
        private static double p_enz_SOC = PlugIn.Parameters.ProportionEnzymeActing;
        private static double mic_to_som = PlugIn.Parameters.FractionMicrobialToSOM;

        // Constants 
        // defaults from Abramoff et al. NOTE:  p[n] indicates the field from the R-code csv input file.
        private static double ea_dep = 64.3413501543551;        //p[1]  #activation energy of SOM depolymerization
        private static double ea_upt = 60.2613470377852;        //p[2]  #activation energy of DOC uptake
        private static double a_dep = 1.0559E+11;         //p[3]   #pre-exponential constant for SOM depolymerization
        private static double a_upt = 1.0794E+11;         //p[4]  #pre-exponential constant for uptake
        private static double km_dep = 0.00246766722424667;        //p[12] #half-saturation constant for SOM depolymerization
        private static double km_upt = 0.289955018042339;        //p[13] #half-saturation constant for DOC uptake
        private static double pconst = 0.481045148509146;        //p[18] #proportion of assimilated C allocated to enzyme production
        private static double qconst = 0.49101177927437;        //p[19] #proportion of assimilated N allocated to enzyme production
        private static double km_o2 = 0.115473482406968;         //p[21] #Michaelis constant for O2
        private static double dgas = 1.63155073368955;          //p[22] #diffusion coefficient for O2 in air
        private static double dliq = 3.13540523153298;          //p[23] #diffusion coefficient for unprotected SOM and DOM in liquid
        private static double o2airfrac = 0.202587071980532;     //p[24] #volume fraction of O2 in air
        private static double saturation = 0.5;     //p[29] #saturation level (was: sat)
        public static double r = 0.008314;                 // gas constant

        // CN Ratios
        public static double CN_DOCN = 27.6;         //p[9] #C:N of soil (was: cns)  // Rob: only used for initializing DON
        private static double cn_microbial = 13;  //p[10] #C:N of microbial biomass, Vogel
        private static double cn_enzymes = 3;    //p[11] #C:N of enzymes (was: cne)

        //private static double frac = 0.000378981;          //p[5]  #fraction of unprotected SOM, Magill et al. 2000
        //private static double frac = 0.00014;          //p[5]  #fraction of unprotected SOM, Vogel
        // p[6] 0.000466501;  NOT USED
        // p[7] 0.000486085;  NOT USED
        //private static double cn_litter = 48.81375915;     //p[8] #C:N of litter (was: cnl) [Not used]
        //private static double cn_microbial = 9.803885526;  //p[10] #C:N of microbial biomass (was: cnm)
        //private static double r_ecloss = 0.000988870450524585;       //p[14] #enzyme turnover rate
        //private static double r_death = 0.000150496216739791;       //p[15] #microbial turnover rate
        //private static double r_death = 0.00150496216739791;       //p[15] #microbial turnover rate
        //private static double c_use_efficiency = 0.31;           //p[16] #carbon use efficiency (was: cue)
        //private static double p_enz_SOC = 0.508531920004333;      //p[17] #proportion of enzyme pool acting on SOC (was: a)
        //private static double mic_to_som = 0.493381458516761;    //p[20] #fraction of dead microbial biomass allocated to SOM
        //private double bulk_density = 0.75743956;         //p[25] #bulk density (was: bd)
        //private double particle_density = 2.50156948;     //p[26] #particle density (was: pd)
        //private static double soilMoistureA = -1.92593874;  //p[27]
        //private static double soilMoistureA = 0.0001;  //p[27]
        //private static double soilMoistureB = 1.0;  //p[28]                     


        public static void Decompose(int Year, int Month, ActiveSite site)
        {

            double month_to_hours =  24.0 * (double) AnnualClimate.DaysInMonth(Month, Year);
            double mg_to_g = 1000;

            double SoilT = SiteVars.SoilTemperature[site];
            double SoilMoisture = SiteVars.SoilWaterContent[site];  //ML: This may be a volume of water while DAMM wants a VSWC
            double bulk_density = SiteVars.SoilBulkDensity[site];
            double particle_density = SiteVars.SoilParticleDensity[site];

            double LitterCinput = (SiteVars.SoilPrimary[site].MonthlyCarbonInputs * 0.5) / month_to_hours * mg_to_g;
            double LitterNinput = SiteVars.SoilPrimary[site].MonthlyNitrogenInputs / month_to_hours * mg_to_g;
            double DOCinput = (SiteVars.SoilPrimary[site].MonthlyCarbonInputs * 0.5) / month_to_hours * mg_to_g;
            double c_loss = 0.0;

            //PlugIn.ModelCore.UI.WriteLine(" SoilLayer:  Month={0}, SoilT={1:0.00}, SoilMoisture={2:0.00}, LitterC={3:0.00}, LitterN={4:0.00}", Month, SoilT, SoilMoisture, LitterCinput, LitterNinput);

            // seasonal DOC input:
            //double A1 = 0.0005; //       #seasonal amplitude
            //double w1 = 2 * Math.PI / 4559;
            //double dref = 0.0005;     //#reference input
            //double t = 1.0;
            // double DOC_input = dref + A1 * Math.Sin(w1 * t - Math.PI / 2);  // # mg C cm-3 soil 

            for (int hours = 0; hours < month_to_hours; hours++)
            {

                double SOC = SiteVars.SoilPrimary[site].Carbon * mg_to_g;
                double SON = SiteVars.SoilPrimary[site].Nitrogen * mg_to_g;
                double DOC = SiteVars.SoilPrimary[site].DOC * mg_to_g;
                double DON = SiteVars.SoilPrimary[site].DON * mg_to_g;
                double microbial_C = SiteVars.SoilPrimary[site].MicrobialCarbon * mg_to_g;
                double microbial_N = SiteVars.SoilPrimary[site].MicrobialNitrogen * mg_to_g;
                double enzymatic_concentration = SiteVars.SoilPrimary[site].EnzymaticConcentration * mg_to_g;

                //PlugIn.ModelCore.UI.WriteLine(" Initial:  Month={0}, SOC={1:0.0}, SON:{2:0.0}, DOC={3:0.0}, DON={4:0.0}, MicC={5:0.0}, MicN={6:0.0}, EC={7:0.0}", Month, SOC, SON, DOC, DON, microbial_C, microbial_N, enzymatic_concentration);
                //double cn_microbial = microbial_C / microbial_N;
                double cn_litter = LitterCinput / LitterNinput;

                double porosity = 1.0 - (bulk_density / particle_density);                                      //calculate porosity                
                                                                                                                //double soilm = -soilMoistureA + soilMoistureB * SoilMoisture;                              //calculate soil moisture scalar, omit gave negative values         
                double soilm = SoilMoisture;                                                                //calculate soil moisture        
                soilm = (soilm > saturation) ? saturation : soilm;                                          //set upper bound on soil moisture (saturation)
                soilm = (soilm < 0.1) ? 0.1 : soilm;                                                        //set lower bound on soil moisture               
                double o2 = dgas * o2airfrac * Math.Pow((porosity - soilm), (4.0 / 3.0));                   //calculate oxygen concentration
                double sol_soc = dliq * Math.Pow(soilm, 3) * frac * SOC;
                double sol_son = dliq * Math.Pow(soilm, 3) * frac * SON;                                    //calculate unprotected SON
                double vmax_dep = a_dep * Math.Exp(-ea_dep / (r * (SoilT + 273)));                          //calculate maximum depolymerization rate                
                double vmax_upt = a_upt * Math.Exp(-ea_upt / (r * (SoilT + 273)));                          //calculate maximum depolymerization rate               
                double upt_c = microbial_C * vmax_upt * DOC / (km_upt + DOC) * o2 / (km_o2 + o2);           //calculate DOC uptake
                double c_mineralization = upt_c * (1.0 - c_use_efficiency);                                   //calculate initial C mineralization               
                double upt_n = microbial_N * vmax_upt * DON / (km_upt + DON) * o2 / (km_o2 + o2);           //calculate DON uptake
                double death_c = r_death * Math.Pow(microbial_C, 2.0);                                        //calculate density-dependent microbial C turnover
                double death_n = r_death * Math.Pow(microbial_N, 2.0);                                        //calculate density-dependent microbial N turnover                
                double enz_c = pconst * c_use_efficiency * upt_c;                                           //calculate potential enzyme C production
                double enz_n = qconst * upt_n;                                                              //calculate potential enzyme N production
                double eprod = (enz_c / cn_enzymes >= enz_n) ? enz_n : (enz_c / cn_enzymes);                //calculate actual enzyme based on Liebig's Law

                double growth_c = (1 - pconst) * (upt_c * c_use_efficiency) + enz_c - cn_enzymes * eprod;   //calculate potential microbial biomass C growth
                double growth_n = (1 - qconst) * upt_n + enz_n - eprod;                                     //calculate potential microbial biomass N growth
                double growth = (growth_c / cn_microbial >= growth_n) ? growth_n : (growth_c / cn_microbial); //calculate actual microbial biomass growth based on Liebig's Law of the minimum (Schimel & Weintraub 2003 SBB)

                double overflow = growth_c - cn_microbial * growth;                         //calculate overflow metabolism of C
                double nmin = growth_n - growth;                                            //calculate N mineralization

                double dmic_c = (cn_microbial * growth) - death_c;                            //calculate change in microbial C pool
                double dmic_n = growth - death_n;                                           //calculate change in microbial N pool

                double eloss = r_ecloss * enzymatic_concentration;                          //calculate enzyme turnover
                double dec = eprod - eloss;                                                 //calculate change in enzyme pool

                double decom_c = vmax_dep * p_enz_SOC * enzymatic_concentration * sol_soc / (km_dep + sol_soc + enzymatic_concentration);     //calculate depolymerization of SOC using ECA kinetics (Tang 2015 GMD)
                double decom_n = vmax_dep * (1 - p_enz_SOC) * enzymatic_concentration * sol_son / (km_dep + sol_son + enzymatic_concentration); //calculate depolymerization of SON using ECA kinetics 

                double dsoc = LitterCinput + (death_c * mic_to_som) - decom_c;                    //calculate change in SOC pool
                double dson = LitterNinput + (death_n * mic_to_som) - decom_n;                    //calculate change in SON pool

                //PlugIn.ModelCore.UI.WriteLine(" DSOC:  Month={0}, LitterInput={1:0.00}, MicrobialDeath={2:0.00}, Decom_C={3:0.00}", Month, LitterCinput, (death_c * mic_to_som), decom_c);
                //PlugIn.ModelCore.UI.WriteLine(" DSON:  Month={0}, LitterInput={1:0.00}, MicrobialDeath={2:0.00}, Decom_N={3:0.00}", Month, LitterNinput, (death_n * mic_to_som), decom_n);

                double ddoc = DOCinput + decom_c + (death_c * (1.0 - mic_to_som)) + (cn_enzymes * eloss) - upt_c; //calculate change in DOC pool
                double ddon = (DOCinput / (DOC / DON)) + decom_n + (death_n * (1.0 - mic_to_som)) + eloss - upt_n; //calculate change in DON pool

                SiteVars.SoilPrimary[site].Carbon += dsoc / mg_to_g;
                SiteVars.SoilPrimary[site].Nitrogen += dson / mg_to_g;
                SiteVars.SoilPrimary[site].DOC = Math.Max(SiteVars.SoilPrimary[site].DOC + (ddoc / mg_to_g), 0.001);
                SiteVars.SoilPrimary[site].DON = Math.Max(SiteVars.SoilPrimary[site].DON + (ddon / mg_to_g), 0.0001);
                SiteVars.SoilPrimary[site].MicrobialCarbon += dmic_c / mg_to_g;
                SiteVars.SoilPrimary[site].MicrobialNitrogen += dmic_n / mg_to_g;
                SiteVars.SoilPrimary[site].EnzymaticConcentration += dec / mg_to_g;

                //double dcout = cmin + overflow;                                             //calculate C efflux

                c_loss += (c_mineralization + overflow) / mg_to_g; //calculate C efflux
            }

            SiteVars.SoilPrimary[site].Respiration(c_loss, site);  // RMS: 11-13-2019: This causes negative mineral N.  

            SiteVars.SoilPrimary[site].MonthlyCarbonInputs = 0.0;  // Done with these now.
            SiteVars.SoilPrimary[site].MonthlyNitrogenInputs = 0.0;  // Done with these now.

            double cLeached = 0.0;  // Carbon leached to a stream

            if (SiteVars.WaterMovement[site] > 0.0)  //Volume of water moving-ML.  Used to be an index of water movement that indicates saturation (amov)
            {

                double leachTextureEffect = OtherData.OMLeachIntercept + OtherData.OMLeachSlope * SiteVars.SoilPercentSand[site];

                double indexWaterMovement = SiteVars.WaterMovement[site] / (SiteVars.SoilDepth[site] * SiteVars.SoilFieldCapacity[site]);

                cLeached = c_loss * leachTextureEffect * indexWaterMovement;

                //Partition and schedule C flows 
                if (cLeached > SiteVars.SoilPrimary[site].Carbon)
                    cLeached = SiteVars.SoilPrimary[site].Carbon;

                //round these to avoid unexpected behavior
                SiteVars.SoilPrimary[site].Carbon = Math.Round((SiteVars.SoilPrimary[site].Carbon - cLeached));
                SiteVars.Stream[site].Carbon = Math.Round((SiteVars.Stream[site].Carbon + cLeached));

                // Compute and schedule N flows and update mineralization accumulators
                // Need to use the ratio for som1 for organic leaching
                double ratioCN_SoilPrimary = SiteVars.SoilPrimary[site].Carbon / SiteVars.SoilPrimary[site].Nitrogen;
                double orgflow = cLeached / ratioCN_SoilPrimary;

                SiteVars.SoilPrimary[site].Nitrogen -= orgflow;
                SiteVars.Stream[site].Nitrogen += orgflow;

                SiteVars.MonthlyStreamN[site][Main.Month] += orgflow;
            }

        }
        public void Respiration(double co2loss, ActiveSite site)
        {
            // Compute flows associated with microbial respiration.

            // Input:
            //  co2loss = CO2 loss associated with decomposition

            // Mineralization associated with respiration is proportional to the N fraction.
            // 12/26/2019: Commented out N dynamics below... creates surge in mineral N.
            //double mineralNFlow = co2loss * this.Nitrogen / this.Carbon;

            //if (mineralNFlow > this.Nitrogen)
            //{
            //    //if((mineralNFlow - this.Nitrogen) > 0.01)
            //    //{
            //    //    PlugIn.ModelCore.UI.WriteLine("RESPIRATION for layer {0} {1}:  Mineral N flow exceeds layer Nitrogen.", this.Name, this.Type);
            //    //    PlugIn.ModelCore.UI.WriteLine("  MineralNFlow={0:0.000}, this.Nitrogen ={0:0.000}", mineralNFlow, this.Nitrogen);
            //    //    PlugIn.ModelCore.UI.WriteLine("  CO2 loss={0:0.000}, this.Carbon={0:0.000}", co2loss, this.Carbon);
            //    //    PlugIn.ModelCore.UI.WriteLine("  Site R/C: {0}/{1}.", site.Location.Row, site.Location.Column);
            //    //}
            //    mineralNFlow = this.Nitrogen;
            //    co2loss = this.Carbon;
            //}

            if (co2loss > this.Carbon)
                co2loss = this.Carbon;

            //round these to avoid unexpected behavior
            SiteVars.SourceSink[site].Carbon = Math.Round((SiteVars.SourceSink[site].Carbon + co2loss));

            //Add lost CO2 to monthly heterotrophic respiration
            SiteVars.MonthlyResp[site][Main.Month] += co2loss;

            //SiteVars.MineralN[site] += mineralNFlow;  // 12/26/2019: Causes surge of Mineral N

            return;
        }
    }
}
