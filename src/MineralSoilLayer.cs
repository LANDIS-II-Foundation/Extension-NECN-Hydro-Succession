//  Author: Robert Scheller, Melissa Lucash

using Landis.Core;
using Landis.SpatialModeling;
using Landis.Utilities;
using System;

namespace Landis.Extension.Succession.NECN
{
    /// <summary>
    /// </summary>
    public class MineralSoilLayer 
    {
        
        public static void Decompose(ActiveSite site)
        {
            
            IEcoregion ecoregion = PlugIn.ModelCore.Ecoregion[site];
            
            double som2c = SiteVars.MineralSoil[site].Carbon;

            if (som2c < 0.0000001)
                return;
            //
            // Determine C/N ratios for flows to SOM1
            //double ratioCNto1 = Layer.BelowgroundDecompositionRatio(site,
            //                        OtherData.MinCNenterSOM1, 
            //                        OtherData.MaxCNenterSOM1,
            //                        OtherData.MinContentN_SOM1);

            double anerb = SiteVars.AnaerobicEffect[site];  

            // Compute total C flow out of SOM2C
            double totalCflow = som2c
                            * SiteVars.DecayFactor[site]
                            * PlugIn.Parameters.DecayRateMineralSoil
                            * anerb //impact of soil anaerobic conditions
                            * OtherData.MonthAdjust;
            ////PlugIn.ModelCore.UI.WriteLine("som2c={0:0.00}, decayFactor={1:0.00}, decayRateSOM2={2:0.00}, anerb={3:0.00}, monthAdj={4:0.00}", som2c, SiteVars.DecayFactor[site], ClimateRegionData.DecayRateSOM2[ecoregion], anerb, OtherData.MonthAdjust);

            //// If SOM2 can decompose to SOM1, it will also go to SOM3.
            //// If it can't go to SOM1, it can't decompose at all.

            //if (SiteVars.SOM2[site].DecomposePossible(ratioCNto1, SiteVars.MineralN[site]))
            //    //PlugIn.ModelCore.UI.WriteLine("DecomposePoss.  MineralN={0:0.00}.", SiteVars.MineralN[site]);
            //{

            //CO2 loss - Compute and schedule respiration flows
            double co2loss = totalCflow * PlugIn.Parameters.FractionMineralSoilToCO2;
            double netCFlow = totalCflow - co2loss;
            SiteVars.MineralSoil[site].Respiration(co2loss, site);
            //PlugIn.ModelCore.UI.WriteLine("AfterTransferto.  MineralN={0:0.00}.", SiteVars.MineralN[site]);

            return;
        }
    }
}
