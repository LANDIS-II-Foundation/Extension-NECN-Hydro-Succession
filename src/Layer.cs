//  Author: Robert Scheller, Melissa Lucash

using Landis.Utilities;
using System;
using System.Threading;
using Landis.Core;
using Landis.SpatialModeling;


namespace Landis.Extension.Succession.NECN
{

    public enum LayerName { Leaf, FineRoot, Wood, CoarseRoot, Metabolic, Structural, Mineral, Other }; 
    public enum LayerType {Surface, Soil, Other} 

    /// <summary>
    /// A Century soil model carbon and nitrogen pool.
    /// </summary>
    public class Layer
    {
        private LayerName name;
        private LayerType type;
        private double carbon;
        private double nitrogen;
        private double decayValue;
        private double fractionLignin;
        //private double netMineralization;
        //private double grossMineralization;


        //---------------------------------------------------------------------
        public Layer(LayerName name, LayerType type)
        {
            this.name = name;
            this.type = type;
            this.carbon = 0.0;
            this.nitrogen = 0.0;

            this.decayValue = 0.0;
            this.fractionLignin = 0.0;

            //this.netMineralization = 0.0;
            //this.grossMineralization = 0.0;

        }
        //---------------------------------------------------------------------
        /// <summary>
        /// Layer Name
        /// </summary>
        public LayerName Name
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
        /// Provides an index to LitterTypeTable
        /// </summary>
        public LayerType Type
        {
            get
            {
                return type;
            }
            set
            {
                type = value;
            }
        }
        //---------------------------------------------------------------------

        /// <summary>
        /// Carbon
        /// </summary>
        public double Carbon
        {
            get
            {
                return carbon;
            }
            set
            {
                carbon = value;
            }
        }
        //---------------------------------------------------------------------

        /// <summary>
        /// Nitrogen
        /// </summary>
        public double Nitrogen
        {
            get
            {
                return nitrogen;
            }
            set
            {
                nitrogen = value;
            }
        }
        //---------------------------------------------------------------------

        /// <summary>
        /// Pool decay rate.
        /// </summary>
        public  double DecayValue
        {
            get
            {
                return decayValue;
            }
            set
            {
                decayValue = value;
            }
        }

        //---------------------------------------------------------------------
        /// <summary>
        /// Pool Carbon:Nitrogen Ratio
        /// </summary>
        public  double FractionLignin
        {
            get
            {
                return fractionLignin;
            }
            set
            {
                fractionLignin = value;
            }
        }

        // --------------------------------------------------
        public void DecomposeStructural(ActiveSite site)
        {
            if (this.Carbon < 0.0000001)
                return;

                double anerb = SiteVars.AnaerobicEffect[site];

                if (this.Type == LayerType.Surface) anerb = 1.0; // No anaerobic effect on surface material

                //Compute total C flow out of structural in layer
                double totalCFlow = System.Math.Min(this.Carbon, OtherData.MaxStructuralC)
                                * SiteVars.DecayFactor[site]
                                * OtherData.LitterParameters[(int)this.Type].DecayRateStrucC
                                * anerb
                                * System.Math.Exp(-1.0 * OtherData.LigninDecayEffect * this.FractionLignin)
                                * OtherData.MonthAdjust;

                //Decompose structural into SOC with CO2 loss.
                this.DecomposeLignin(totalCFlow, site);
        }
        // --------------------------------------------------
        // Decomposition of compartment lignin
        public void DecomposeLignin(double totalCFlow, ActiveSite site)
        {
            double carbonToOHorizon;    //Net C flow to OHorizon
            double carbonToMineralSoil;    //Net C flow to MineralSoil
            double litterC = this.Carbon; 
            double ratioCN = this.Carbon / this.Nitrogen;


            //See if Layer can decompose to SOM1.
            //If it can decompose to SOM1, it will also go to SOM2.
            //If it can't decompose to SOM1, it can't decompose at all.

            //If Wood Object can decompose
            if (this.DecomposePossible(ratioCN, SiteVars.MineralN[site]))
            {

                // Decompose Wood Object to MineralSoil
                // -----------------------
                // Gross C flow to MineralSoil
                carbonToMineralSoil = totalCFlow * this.FractionLignin;

                //MicrobialRespiration associated with decomposition to som2
                double co2loss = carbonToMineralSoil * OtherData.LigninRespirationRate;

                this.Respiration(co2loss, site);

                //Net C flow to MineralSoil
                carbonToMineralSoil -= co2loss;

                // Partition and schedule C flows 
                if (carbonToMineralSoil > this.Carbon)
                    carbonToMineralSoil = this.Carbon;
                
                SiteVars.MineralSoil[site].Carbon += Math.Round(carbonToMineralSoil, 2);  //rounding to avoid unexpected behavior
                SiteVars.MineralSoil[site].Nitrogen += this.TransferNitrogen(carbonToMineralSoil, litterC, ratioCN, site);

                // ----------------------------------------------
                // Decompose Wood Object to OHorizon
                // Gross C flow to OHorizon

                carbonToOHorizon = totalCFlow - carbonToMineralSoil;

                //MicrobialRespiration associated with decomposition to OHorizon
                if (this.Type == LayerType.Surface)
                    co2loss = carbonToOHorizon * OtherData.StructuralToCO2Surface;
                else
                    co2loss = carbonToOHorizon * OtherData.StructuralToCO2Soil;

                this.Respiration(co2loss, site);

                //Net C flow to SOM1
                carbonToOHorizon -= co2loss;

                SiteVars.OHorizon[site].MonthlyCarbonInputs += Math.Round(carbonToOHorizon, 2);  //rounding to avoid unexpected behavior
                SiteVars.OHorizon[site].MonthlyNitrogenInputs += this.TransferNitrogen(carbonToOHorizon, litterC, ratioCN, site);


                //if(this.Type == LayerType.Surface)
                //{
                //    this.TransferCarbon(SiteVars.SOM1surface[site], carbonToSOM1);
                //    this.TransferNitrogen(SiteVars.SOM1surface[site], carbonToSOM1, litterC, ratioCN, site);
                //}
                //else
                //{
                //this.TransferCarbon(SiteVars.OHorizon[site], carbonToOHorizon);
                //    this.TransferNitrogen(SiteVars.OHorizon[site], carbonToOHorizon, litterC, ratioCN, site);
                //}
            }
            //PlugIn.ModelCore.UI.WriteLine("Decompose2.  MineralN={0:0.00}.", SiteVars.MineralN[site]);
            return;
        }

        //---------------------------------------------------------------------
        public void DecomposeMetabolic(ActiveSite site)
        {
            double litterC = this.Carbon;
            double anerb = SiteVars.AnaerobicEffect[site];
            
            // RMS:  Assume that OHorizon is 10cm, the rest is MineralSoil.  
            double fractionOHorizon = Math.Max(10.0 / SiteVars.SoilDepth[site], 1.0);
            double fractionMineralSoil = 1.0 - fractionOHorizon;

            if (litterC > 0.0000001)
                return;

            // Determine C/N ratios for flows to OHorizon
            double ratioCNtoSoils = 0.0;
            double co2loss = 0.0;

            // Compute ratios for surface  metabolic residue
            if (this.Type == LayerType.Surface)
                ratioCNtoSoils = Layer.AbovegroundDecompositionRatio(this.Nitrogen, litterC);

            //Compute ratios for soil metabolic residue
            else
                ratioCNtoSoils = Layer.BelowgroundDecompositionRatio(site,
                                    OtherData.MinCNenterOHorizon,
                                    OtherData.MaxCNenterOHorizon,
                                    OtherData.MinContentN_OHorizon);

            //Compute total C flow out of metabolic layer
            double totalCFlow = litterC
                            * SiteVars.DecayFactor[site]
                            * OtherData.LitterParameters[(int)this.Type].DecayRateMetabolicC
                            * OtherData.MonthAdjust;

            //Effect of soil anerobic conditions: 
            if (this.Type == LayerType.Soil) totalCFlow *= anerb;

            // Add DOC from litter to OHorizon
            double flowToDOC = totalCFlow * PlugIn.Parameters.FractionLitterDecayToDOC * fractionOHorizon;
            double flowToDON = flowToDOC / ratioCNtoSoils;

            this.Carbon -= Math.Min(flowToDOC, this.Carbon);
            SiteVars.OHorizon[site].DOC += flowToDOC;
            this.Nitrogen -= Math.Min(flowToDON, this.Nitrogen);
            SiteVars.OHorizon[site].DON += flowToDON;

            //Make sure metabolic C does not go negative.
            if (totalCFlow > litterC)
                totalCFlow = litterC;

            //If decomposition can occur,
            if (this.DecomposePossible(ratioCNtoSoils, SiteVars.MineralN[site]))
            {
                //CO2 loss
                if (this.Type == LayerType.Surface)
                    co2loss = totalCFlow * OtherData.MetabolicToCO2Surface;
                else
                    co2loss = totalCFlow * OtherData.MetabolicToCO2Soil;

                this.Respiration(co2loss, site);


                //Decompose metabolic into SOC / SON
                double carbonToSoils = totalCFlow - co2loss;

                if (carbonToSoils > litterC && PlugIn.Verbose)
                    PlugIn.ModelCore.UI.WriteLine("   ERROR:  Decompose Metabolic:  netCFlow={0:0.000} > layer.Carbon={0:0.000}.", carbonToSoils, this.Carbon);

                SiteVars.OHorizon[site].MonthlyCarbonInputs += Math.Round(carbonToSoils * fractionOHorizon, 2);  //rounding to avoid unexpected behavior
                SiteVars.OHorizon[site].MonthlyNitrogenInputs += this.TransferNitrogen(carbonToSoils * fractionOHorizon, litterC, ratioCNtoSoils, site);

                // RMS:  Some portion goes directly to MineralSoil, depending on source and depth relative to rooting dept.
                SiteVars.MineralSoil[site].Carbon += Math.Round(carbonToSoils * fractionMineralSoil, 2);  //rounding to avoid unexpected behavior
                SiteVars.MineralSoil[site].Nitrogen += this.TransferNitrogen(carbonToSoils * fractionMineralSoil, litterC, ratioCNtoSoils, site);
                
                //this.TransferCarbon(SiteVars.OHorizon[site], netCFlow);
                //this.TransferNitrogen(SiteVars.OHorizon[site], netCFlow, litterC, ratioCNtoOHorizon, site);

                // -- CARBON AND NITROGEN ---------------------------
                // Partition and schedule C flows
                // Compute and schedule N flows and update mineralization accumulators.
                //if((int) this.Type == (int) LayerType.Surface)
                //{
                //    this.TransferCarbon(SiteVars.SOM1surface[site], netCFlow);
                //    this.TransferNitrogen(SiteVars.SOM1surface[site], netCFlow, litterC, ratioCNtoSOM1, site);
                //    //PlugIn.ModelCore.UI.WriteLine("DecomposeMetabolic.  MineralN={0:0.00}.", SiteVars.MineralN[site]);
                //}
                //else
                //{
                //}

            }
        }
        ////---------------------------------------------------------------------
        //public void TransferCarbon(Layer destination, double netCFlow)
        //{
        //    if (netCFlow < 0)
        //    {
        //        //PlugIn.ModelCore.UI.WriteLine("NEGATIVE C FLOW!  Source: {0},{1}; Destination: {2},{3}.", this.Name, this.Type, destination.Name, destination.Type);
        //    }

        //    if (netCFlow > this.Carbon)
        //        netCFlow = this.Carbon;
        //    //PlugIn.ModelCore.UI.WriteLine("C FLOW EXCEEDS SOURCE!  Source: {0},{1}; Destination: {2},{3}.", this.Name, this.Type, destination.Name, destination.Type);

        //    //round these to avoid unexpected behavior
        //    this.Carbon = Math.Round((this.Carbon - netCFlow), 2);
        //    destination.MonthlyCarbonInputs += Math.Round(netCFlow, 2);
        //}

        public double TransferNitrogen(double CFlow, double totalC, double ratioCNtoDestination, ActiveSite site)
        {
            double mineralNFlow = 0.0;

            //...N flow is proportional to C flow.
            double NFlow = this.Nitrogen * CFlow / totalC;

            //...This was added to avoid a 0/0 error on the pc.
            if (CFlow <= 0.0 || NFlow <= 0.0)
            {
                return 0.0;
            }

            if ((NFlow - this.Nitrogen) > 0.01 && PlugIn.Verbose)
            {
                PlugIn.ModelCore.UI.WriteLine("  Transfer N:  N flow > source N.");
                PlugIn.ModelCore.UI.WriteLine("     NFlow={0:0.000}, SourceN={1:0.000}", NFlow, this.Nitrogen);
                PlugIn.ModelCore.UI.WriteLine("     CFlow={0:0.000}, totalC={1:0.000}", CFlow, totalC);
                PlugIn.ModelCore.UI.WriteLine("     this.Name={0}, this.Type={1}", this.Name, this.Type);
                //PlugIn.ModelCore.UI.WriteLine("     dest.Name  ={0}, dest.Type  ={1}", destination.Name, destination.Type);
                PlugIn.ModelCore.UI.WriteLine("     ratio CN to dest={0}", ratioCNtoDestination);
            }

            double nitrogenInputs = 0.0;

            //...If C/N of Box A > C/N of new material entering Box B, IMMOBILIZATION occurs.
            //   Compute the amount of N immobilized.
            //   since  ratioCNtoDestination = netCFlow / (Nflow + immobileN),
            //   where immobileN is the extra N needed from the mineral pool
            if ((CFlow / NFlow) > ratioCNtoDestination)
            {
                double immobileN = (CFlow / ratioCNtoDestination) - NFlow;
                
                //...Schedule flow from Box A to Box B
                this.Nitrogen -= NFlow;
                //destination.Nitrogen += NFlow;
                //destination.MonthlyNitrogenInputs += NFlow;
                nitrogenInputs = NFlow;

                if(PlugIn.Verbose)
                    PlugIn.ModelCore.UI.WriteLine("TransferNitrogen: Before immobilization: MineralN={0:0.00}, ImmobileN={1:0.000}.", SiteVars.MineralN[site],immobileN);

                // Schedule flow from mineral pool to immobileN
                // Don't allow mineral N to go to zero or negative.- ML

                if (immobileN > SiteVars.MineralN[site])
                    immobileN = SiteVars.MineralN[site] - 0.01; //leave some small amount of mineral N


                SiteVars.MineralN[site] -= immobileN;
                if (PlugIn.Verbose)
                    PlugIn.ModelCore.UI.WriteLine("TransferNitrogen: After immobilization: MineralN={0:0.00}.", SiteVars.MineralN[site]);
                
                //destination.Nitrogen += immobileN;
                nitrogenInputs += immobileN;

                //PlugIn.ModelCore.UI.WriteLine("AdjustImmobil.  MineralN={0:0.00}.", SiteVars.MineralN[site]);
                //PlugIn.ModelCore.UI.WriteLine("   TransferN immobileN={0:0.000}, C={1:0.000}, N={2:0.000}, ratioCN={3:0.000}.", immobileN, CFlow, NFlow, ratioCNtoDestination);
                //PlugIn.ModelCore.UI.WriteLine("     source={0}-{1}, destination={2}-{3}.", this.Name, this.Type, destination.Name, destination.Type);

                //...Return mineralization value.
                mineralNFlow = -1 * immobileN;
                //PlugIn.ModelCore.UI.WriteLine("MineralNflow.  MineralN={0:0.00}.", SiteVars.MineralN[site]);
            }
            else

            // MINERALIZATION occurs 
            {
                double mineralizedN = (CFlow / ratioCNtoDestination);

                this.Nitrogen -= mineralizedN;
                //destination.Nitrogen += mineralizedN;
                nitrogenInputs += mineralizedN;

                //...Schedule flow from Box A to mineral pool

                mineralNFlow = NFlow - mineralizedN;

                //if ((mineralNFlow - this.Nitrogen) > 0.01 && PlugIn.Verbose) 
                if(PlugIn.Verbose)
                {
                    PlugIn.ModelCore.UI.WriteLine("  Nitrogen Mineralization");
                    PlugIn.ModelCore.UI.WriteLine("  Transfer N mineralization:  mineralN > source N.");
                    PlugIn.ModelCore.UI.WriteLine("     MineralNFlow={0:0.000}, SourceN={1:0.000}", mineralNFlow, this.Nitrogen);
                    PlugIn.ModelCore.UI.WriteLine("     CFlow={0:0.000}, totalC={1:0.000}", CFlow, totalC);
                    PlugIn.ModelCore.UI.WriteLine("     this.Name={0}, this.Type={1}", this.Name, this.Type);
                    //PlugIn.ModelCore.UI.WriteLine("     dest.Name  ={0}, dest.Type  ={1}", destination.Name, destination.Type);
                    PlugIn.ModelCore.UI.WriteLine("     ratio CN to dest={0}", ratioCNtoDestination);
                }

                //this.Nitrogen -= mineralNFlow;

                SiteVars.MineralN[site] += mineralNFlow;
                if (mineralNFlow > 3.0 && PlugIn.Verbose)
                    PlugIn.ModelCore.UI.WriteLine("Layer.TransferNitrogen: N Mineralization = {0:0.00})", mineralNFlow);

            }

            if (mineralNFlow > 0)
                SiteVars.GrossMineralization[site] += mineralNFlow;

            return nitrogenInputs;
        }

        public void Respiration(double co2loss, ActiveSite site)
        {
            // Compute flows associated with microbial respiration.
            //  co2loss = CO2 loss associated with decomposition

            //  Output:
            //  carbonSourceSink = C source/sink
            //  grossMineralization = gross mineralization
            //  netMineralization = net mineralization for layer N

            //c...Mineralization associated with respiration is proportional to the N fraction.
            if (this.Nitrogen < 0.000001 || this.Carbon < 0.00001)  // Not enough C or N to respire
                return;

            double mineralNFlow = co2loss * this.Nitrogen / this.Carbon;
            if (PlugIn.Verbose)
            {
                PlugIn.ModelCore.UI.WriteLine("Layer.Respiration: Source:  this.Name={0}, this.Type={1}", this.Name, this.Type);
                PlugIn.ModelCore.UI.WriteLine("Layer.Respiration: this.Nitrogen= {0:0.000}, this.Carbon={1:0.00}", this.Nitrogen, this.Carbon);
            }

            if (mineralNFlow > this.Nitrogen)
            {
                if (PlugIn.Verbose)
                {
                    PlugIn.ModelCore.UI.WriteLine("RESPIRATION for layer {0} {1}:  Mineral N flow exceeds layer Nitrogen.", this.Name, this.Type);
                    PlugIn.ModelCore.UI.WriteLine("  MineralNFlow={0:0.000}, this.Nitrogen ={0:0.000}", mineralNFlow, this.Nitrogen);
                    PlugIn.ModelCore.UI.WriteLine("  CO2 loss={0:0.000}, this.Carbon={0:0.000}", co2loss, this.Carbon);
                    PlugIn.ModelCore.UI.WriteLine("  Site R/C: {0}/{1}.", site.Location.Row, site.Location.Column);
                }
                mineralNFlow = this.Nitrogen;
                co2loss = this.Carbon;
            }

            if (co2loss > this.Carbon)
                co2loss = this.Carbon;

            //round these to avoid unexpected behavior
            this.Carbon = Math.Round((this.Carbon - co2loss));
            SiteVars.SourceSink[site].Carbon = Math.Round((SiteVars.SourceSink[site].Carbon + co2loss));

            //Add loss CO2 to monthly heterotrophic respiration
            SiteVars.MonthlyResp[site][Main.Month] += co2loss;

            this.Nitrogen -= mineralNFlow;
            SiteVars.MineralN[site] += mineralNFlow;
            if (mineralNFlow > 3.0 && PlugIn.Verbose)
                PlugIn.ModelCore.UI.WriteLine("Layer.Respiration: N Mineralization = {0:0.00})", mineralNFlow);

            if (PlugIn.Verbose)
            {
                PlugIn.ModelCore.UI.WriteLine("Layer.Respiration: Source:  this.Name={0}, this.Type={1}", this.Name, this.Type);
                PlugIn.ModelCore.UI.WriteLine("Layer.Respiration: MineralNflow= {0:0.000}, co2loss={1:0.00}", mineralNFlow, co2loss);
            }
           
            //c...Update gross mineralization
            // this.GrossMineralization += mineralNFlow;
            if (mineralNFlow > 0)
                SiteVars.GrossMineralization[site] += mineralNFlow;

            //c...Update net mineralization
            //this.NetMineralization += mineralNFlow;

            return;
        }

        public bool DecomposePossible(double ratioCNnew, double mineralN)
        {

            //c...Determine if decomposition can occur.

            bool canDecompose = true;

            //c...If there is no available mineral N
            if (mineralN < 0.0000001)
            {

                // Compare the C/N of new material to the C/N of the layer if C/N of
                // the layer > C/N of new material
                if (this.Carbon / this.Nitrogen > ratioCNnew)
                {

                    // Immobilization is necessary and the stuff in Box A can't
                    // decompose to Box B.
                    canDecompose = false;
                }
            }

            // If there is some available mineral N, decomposition can
            // proceed even if mineral N is driven negative in
            // the next time step.

            return canDecompose;

        }

        public void AdjustLignin(double inputC, double inputFracLignin)
        {
            //c...Adjust the fraction of lignin in structural C when new material
            //c...  is added.

            //c    oldc  = grams C in structural before new material is added
            //c    frnew = fraction of lignin in new structural material
            //c    addc  = grams structural C being added

            //c    fractl comes in as fraction of lignin in structural before new
            //c           material is added; goes out as fraction of lignin in
            //c           structural with old and new combined.

            //c...oldlig  = grams of lignin in existing residue
            double oldlig = this.FractionLignin * this.Carbon;//totalC;

            //c...newlig = grams of lignin in new residue
            double newlig = inputFracLignin * inputC;

            //c...Compute lignin fraction in combined residue
            double newFraction = (oldlig + newlig) / (this.Carbon + inputC);

            this.FractionLignin = newFraction;

            return;
        }

        public void AdjustDecayRate(double inputC, double inputDecayRate)
        {
            //c...oldlig  = grams of lignin in existing residue
            double oldDecayRate = this.DecayValue * this.Carbon;

            //c...newlig = grams of lignin in new residue
            double newDecayRate = inputDecayRate * inputC;

            //c...Compute decay rate in combined residue
            this.DecayValue = (oldDecayRate + newDecayRate) / (inputC + this.Carbon);

            return;
        }


        public static double BelowgroundDecompositionRatio(ActiveSite site, double minCNenter, double maxCNenter, double minContentN)
        {
            //BelowGround Decomposition RATio computation.
            double bgdrat = 0.0;

            //Determine ratio of C/N of new material entering 'Box B'.
            //Ratio depends on available N

            double mineralN = SiteVars.MineralN[site];

            if (mineralN <= 0.0)
                bgdrat = maxCNenter;  // Set ratio to maximum allowed (HIGHEST carbon, LOWEST nitrogen)
            else if (mineralN > minContentN)
                bgdrat = minCNenter;  //Set ratio to minimum allowed
            else
                bgdrat = (1.0 - (mineralN / minContentN)) * (maxCNenter - minCNenter)
                    + minCNenter;

            return bgdrat;
        }

        public static double AbovegroundDecompositionRatio(double abovegroundN, double abovegroundC)
        {       

            double Ncontent, agdrat;
            double biomassConversion = 2.0;
            
            // cemicb = slope of the regression line for C/N of som1
            double cemicb = (OtherData.MinCNSurfMicrobes - OtherData.MaxCNSurfMicrobes) / OtherData.MinNContentCNSurfMicrobes;


            //The C/E ratios for structural and wood can be computed once;
            //they then remain fixed throughout the run.  The ratios for
            //metabolic and som1 may vary and must be recomputed each time step

            if ((abovegroundC * biomassConversion) <= 0.00000000001)  Ncontent = 0.0;
            else  Ncontent = abovegroundN / (abovegroundC * biomassConversion);

            //tca is multiplied by biomassConversion to give biomass

            if (Ncontent > OtherData.MinNContentCNSurfMicrobes)
                agdrat = OtherData.MinCNSurfMicrobes;
            else
                agdrat = OtherData.MaxCNSurfMicrobes + Ncontent * cemicb;

            return agdrat;
        }

        //---------------------------------------------------------------------
        /// <summary>
        /// Reduces the pool's biomass by a specified percentage.
        /// </summary>
        public void ReduceMass(double percentageLost)
        {
            if (percentageLost < 0.0 || percentageLost > 1.0)
                throw new ArgumentException("Percentage must be between 0% and 100%");

            this.Carbon   = this.Carbon * (1.0 - percentageLost);
            this.Nitrogen   = this.Nitrogen * (1.0 - percentageLost);

            return;
        }

    }
}
