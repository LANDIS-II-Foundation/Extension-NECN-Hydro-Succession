﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Landis.Library.Metadata;

namespace Landis.Extension.Succession.NECN
{
    public class PrimaryLog
    {
            //log.WriteLine("");

        [DataFieldAttribute(Unit = FieldUnits.Year, Desc = "Simulation Year")]
        public int Time {set; get;}

        [DataFieldAttribute(Desc = "Name of Climate Region")]
        public string ClimateRegionName { set; get; }

        [DataFieldAttribute(Desc = "Climate Region Index")]
        public int ClimateRegionIndex { set; get; }

        [DataFieldAttribute(Unit = FieldUnits.Count, Desc = "Number of Sites")]
        public int NumSites { set; get; }

        [DataFieldAttribute(Unit = FieldUnits.g_C_m2, Desc = "Net Ecosystem Exchange C", Format = "0.0")]
        public double NEEC {get; set;}

        [DataFieldAttribute(Unit = FieldUnits.g_C_m2, Desc = "Total Soil Organic Carbon", Format = "0.0")]
        public double SOMTC { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_B_m2, Desc = "Aboveground Biomass", Format = "0.0")]
        public double AGB { get; set; }
        
        //log.Write("AG_NPPC, BG_NPPC, LitterfallC, AgeMortality, ");
        [DataFieldAttribute(Unit = FieldUnits.g_C_m2_yr1, Desc = "Aboveground NPP C", Format = "0.0")]
        public double AG_NPPC { get; set; }
        
        [DataFieldAttribute(Unit = FieldUnits.g_C_m2_yr1, Desc = "Below ground NPP C", Format = "0.0")]
        public double BG_NPPC { get; set; }
        
        [DataFieldAttribute(Unit = FieldUnits.g_C_m2_yr1, Desc = "Litterfall C", Format = "0.0")]
        public double Litterfall { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_B_m2_yr1, Desc = "Age Mortality Biomass", Format = "0.0")]
        public double AgeMortality { get; set; }
        
        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "Mineral N", Format = "0.00")]
        public double MineralN { get; set; }
        
        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "Total N", Format = "0.0")]
        public double TotalN { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_C_m2, Desc = "OHorizon SOC", Format = "0.0")]
        public double C_OHorizon { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_C_m2, Desc = "MineralSoil C", Format = "0.0")]
        public double C_MineralSoil { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_C_m2, Desc = "Live Leaf C", Format = "0.0")]
        public double C_LiveLeaf { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_C_m2, Desc = "Live Fine Root C", Format = "0.0")]
        public double C_LiveFRoot { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_C_m2, Desc = "Live Wood C", Format = "0.0")]
        public double C_LiveWood { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_C_m2, Desc = "Live Coarse Root C", Format = "0.0")]
        public double C_LiveCRoot { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_C_m2, Desc = "Dead Wood C", Format = "0.0")]
        public double C_DeadWood { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_C_m2, Desc = "Dead Coarse Root C", Format = "0.0")]
        public double C_DeadCRoot { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_C_m2, Desc = "Dead Leaf C", Format = "0.0")]
        public double C_DeadLeaf { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_C_m2, Desc = "Dead Fine Root C", Format = "0.0")]
        public double C_DeadFRoot { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "OHorizon SON", Format = "0.0")]
        public double SON_OHorizon { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "OHorizon DON", Format = "0.0")]
        public double DON_OHorizon { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "OHorizon MicN", Format = "0.0")]
        public double MicrobialN_OHorizon { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "MineralSoil N", Format = "0.0")]
        public double N_MineralSoil { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "Leaf N", Format = "0.0")]
        public double N_Leaf { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "Fine Root N", Format = "0.0")]
        public double N_FRoot { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "Wood N", Format = "0.0")]
        public double N_Wood { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "Coarse Root N", Format = "0.0")]
        public double N_CRoot { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "Dead Wood N", Format = "0.0")]
        public double N_DeadWood { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "Dead Coarse Root N", Format = "0.0")]
        public double N_DeadCRoot { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "Dead Leaf N", Format = "0.0")]
        public double N_DeadLeaf { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "Dead Fine Root N", Format = "0.0")]
        public double N_DeadFRoot { get; set; }


        //[DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "Surface Structural Net Mineralization", Format = "0.0")]
        //public double SurfStrucNetMin { get; set; }

        //[DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "Surface Metabolic Net Mineralization", Format = "0.0")]
        //public double SurfMetaNetMin { get; set; }

        //[DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "Soil Structural Net Mineralization", Format = "0.0")]
        //public double SoilStrucNetMin { get; set; }

        //[DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "Soil Metabolic Net Mineralization", Format = "0.0")]
        //public double SoilMetaNetMin { get; set; }

        //log.Write("SOM1surfNetMin, SoilPrimaryNetMin, SOM2NetMin, SOM3NetMin, ");
        //[DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "SOM1 Surface Net Mineralization", Format = "0.0")]
        //public double SOM1surfNetMin { get; set; }

        //[DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "SOM1 Soil Net Mineralization", Format = "0.0")]
        //public double SoilPrimaryNetMin { get; set; }

        //[DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "SOM2 Net Mineralization", Format = "0.0")]
        //public double SOM2NetMin { get; set; }

        //[DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "SOM3 Net Mineralization", Format = "0.0")]
        //public double SOM3NetMin { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "Total Nitrogen Deposition per Timestep", Format = "0.00")]
        public double TotalNdep { get; set; }

        //log.Write("StreamC, StreamN, FireCEfflux, FireNEfflux, ");
        [DataFieldAttribute(Unit = FieldUnits.g_C_m2, Desc = "Stream C", Format = "0.00")]
        public double StreamC { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "Stream N", Format = "0.00")]
        public double StreamN { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_C_m2, Desc = "Fire C Efflux", Format = "0.0")]
        public double FireCEfflux { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "Fire N Efflux", Format = "0.0")]
        public double FireNEfflux { get; set; }
        
        //log.Write("Nuptake, Nresorbed, TotalSoilN, Nvol, avgfrassC,");
        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "N Uptake", Format = "0.0")]
        public double Nuptake { get; set; }
        
        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "N Resorbed", Format = "0.0")]
        public double Nresorbed { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "Total Soil N", Format = "0.0")]
        public double TotalSoilN { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "N Volatilized", Format = "0.00")]
        public double Nvol { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_C_m2, Desc = "Frass C", Format = "0.0")]
        public double FrassC { get; set; }

    }
}
