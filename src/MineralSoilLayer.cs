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

            double anerb = SiteVars.AnaerobicEffect[site];  

            // Compute total C flow out of SOM2C
            double totalCflow = som2c
                            * SiteVars.DecayFactor[site]
                            * PlugIn.Parameters.DecayRateMineralSoil
                            * anerb //impact of soil anaerobic conditions
                            * OtherData.MonthAdjust;
            ////PlugIn.ModelCore.UI.WriteLine("som2c={0:0.00}, decayFactor={1:0.00}, decayRateSOM2={2:0.00}, anerb={3:0.00}, monthAdj={4:0.00}", som2c, SiteVars.DecayFactor[site], ClimateRegionData.DecayRateSOM2[ecoregion], anerb, OtherData.MonthAdjust);

            //CO2 loss - Compute and schedule respiration flows
            double co2loss = totalCflow * PlugIn.Parameters.FractionMineralSoilToCO2;
            double netCFlow = totalCflow - co2loss;
            SiteVars.MineralSoil[site].Respiration(co2loss, site);
            //PlugIn.ModelCore.UI.WriteLine("AfterTransferto.  MineralN={0:0.00}.", SiteVars.MineralSoil[site].Nitrogen);

            // Leaching next -----------------------------------------------------------

            double cLeached = 0.0;  /// Carbon leached to a stream

            if (SiteVars.WaterMovement[site] > 0.0)  //Volume of water moving-ML.  
            {

                double leachTextureEffect = OtherData.OMLeachIntercept + OtherData.OMLeachSlope * SiteVars.SoilPercentSand[site];

                double indexWaterMovement = SiteVars.WaterMovement[site] / (SiteVars.SoilDepth[site] * SiteVars.SoilFieldCapacity[site]);

                cLeached = netCFlow * leachTextureEffect * indexWaterMovement;

                //Partition and schedule C flows 
                if (cLeached > SiteVars.MineralSoil[site].Carbon)
                    cLeached = SiteVars.MineralSoil[site].Carbon;

                //round these to avoid unexpected behavior
                SiteVars.MineralSoil[site].Carbon = Math.Round((SiteVars.MineralSoil[site].Carbon - cLeached));
                SiteVars.Stream[site].Carbon = Math.Round((SiteVars.Stream[site].Carbon + cLeached));

                // Compute and schedule N flows and update mineralization accumulators
                double ratioCN_MineralSoil = SiteVars.MineralSoil[site].Carbon / SiteVars.MineralSoil[site].Nitrogen;
                double orgflow = cLeached / ratioCN_MineralSoil;

                SiteVars.MineralSoil[site].Nitrogen -= orgflow;
                SiteVars.Stream[site].Nitrogen += orgflow;

                SiteVars.MonthlyStreamN[site][Main.Month] += orgflow;
            }


            return;
        }
    }
}
