//  Author: Robert Scheller

using Landis.Core;
using Landis.SpatialModeling;
using Landis.Library.LeafBiomassCohorts;
using Landis.Utilities;
using System.IO;

namespace Landis.Extension.Succession.NECN
{
    /// <summary>
    /// Utility methods.
    /// </summary>
    public static class Util
    {

        /// <summary>
        /// Converts a table indexed by species and ecoregion into a
        /// 2-dimensional array.
        /// </summary>
        //public static T[,] ToArray<T>(Species.AuxParm<Ecoregions.AuxParm<T>> table)
        //{
        //    T[,] array = new T[PlugIn.ModelCore.Ecoregions.Count, PlugIn.ModelCore.Species.Count];
        //    foreach (ISpecies species in PlugIn.ModelCore.Species) {
        //        foreach (IEcoregion ecoregion in PlugIn.ModelCore.Ecoregions) {
        //            array[ecoregion.Index, species.Index] = table[species][ecoregion];
        //        }
        //    }
        //    return array;
        //}
        //---------------------------------------------------------------------

        public static void ReadSoilDepthMap(string path)
        {
            IInputRaster<DoublePixel> map = MakeDoubleMap(path);

 
            using (map)
            {
                DoublePixel pixel = map.BufferPixel;
                foreach (Site site in PlugIn.ModelCore.Landscape.AllSites)
                {
                    map.ReadBufferPixel();
                    int mapValue = (int) pixel.MapCode.Value;
                    if (site.IsActive)
                    {
                            if (mapValue < 0 || mapValue > 300)
                                throw new InputValueException(mapValue.ToString(),
                                                              "Soil depth value {0} is not between {1:0.0} and {2:0.0}. Site_Row={3:0}, Site_Column={4:0}",
                                                              mapValue, 0, 300, site.Location.Row, site.Location.Column);
                        SiteVars.SoilDepth[site] = mapValue;
                    }
                }
            }
        }
        //---------------------------------------------------------------------

        public static void ReadSoilDrainMap(string path)
        {
            IInputRaster<DoublePixel> map = MakeDoubleMap(path);

            using (map)
            {
                DoublePixel pixel = map.BufferPixel;
                foreach (Site site in PlugIn.ModelCore.Landscape.AllSites)
                {
                    map.ReadBufferPixel();
                    double mapValue = pixel.MapCode.Value;
                    if (site.IsActive)
                    {
                        if (mapValue < 0.0 || mapValue > 1.0)
                            throw new InputValueException(mapValue.ToString(),
                                                          "SOil drainage value {0} is not between {1:0.0} and {2:0.0}. Site_Row={3:0}, Site_Column={4:0}",
                                                          mapValue, 0.0, 1.0, site.Location.Row, site.Location.Column);
                        SiteVars.SoilDrain[site] = mapValue;
                    }
                }
            }
        }
        //---------------------------------------------------------------------

        public static void ReadSoilBaseFlowMap(string path)
        {
            IInputRaster<DoublePixel> map = MakeDoubleMap(path);

            using (map)
            {
                DoublePixel pixel = map.BufferPixel;
                foreach (Site site in PlugIn.ModelCore.Landscape.AllSites)
                {
                    map.ReadBufferPixel();
                    double mapValue = pixel.MapCode.Value;
                    if (site.IsActive)
                    {
                        if (mapValue < 0.0 || mapValue > 1.0)
                            throw new InputValueException(mapValue.ToString(),
                                                          "Soil base flow value {0} is not between {1:0.0} and {2:0.0}. Site_Row={3:0}, Site_Column={4:0}",
                                                          mapValue, 0.0, 1.0, site.Location.Row, site.Location.Column);
                        SiteVars.SoilBaseFlowFraction[site] = mapValue;
                    }
                }
            }
        }
        //---------------------------------------------------------------------

        public static void ReadSoilStormFlowMap(string path)
        {
            IInputRaster<DoublePixel> map = MakeDoubleMap(path);

            //try
            //{
            //    map = PlugIn.ModelCore.OpenRaster<DoublePixel>(path);
            //}
            //catch (FileNotFoundException)
            //{
            //    string mesg = string.Format("Error: The file {0} does not exist", path);
            //    throw new System.ApplicationException(mesg);
            //}

            //if (map.Dimensions != PlugIn.ModelCore.Landscape.Dimensions)
            //{
            //    string mesg = string.Format("Error: The input map {0} does not have the same dimension (row, column) as the ecoregions map", path);
            //    throw new System.ApplicationException(mesg);
            //}

            using (map)
            {
                DoublePixel pixel = map.BufferPixel;
                foreach (Site site in PlugIn.ModelCore.Landscape.AllSites)
                {
                    map.ReadBufferPixel();
                    double mapValue = pixel.MapCode.Value;
                    if (site.IsActive)
                    {
                        if (mapValue < 0.0 || mapValue > 1.0)
                            throw new InputValueException(mapValue.ToString(),
                                                          "Soil storm flow value {0} is not between {1:0.0} and {2:0.0}. Site_Row={3:0}, Site_Column={4:0}",
                                                          mapValue, 0.0, 1.0, site.Location.Row, site.Location.Column);
                        SiteVars.SoilStormFlowFraction[site] = mapValue;
                    }
                }
            }
        }
        //---------------------------------------------------------------------

        public static void ReadFieldCapacityMap(string path)
        {
            IInputRaster<DoublePixel> map = MakeDoubleMap(path); 

            //try
            //{
            //    map = PlugIn.ModelCore.OpenRaster<DoublePixel>(path);
            //}
            //catch (FileNotFoundException)
            //{
            //    string mesg = string.Format("Error: The file {0} does not exist", path);
            //    throw new System.ApplicationException(mesg);
            //}

            //if (map.Dimensions != PlugIn.ModelCore.Landscape.Dimensions)
            //{
            //    string mesg = string.Format("Error: The input map {0} does not have the same dimension (row, column) as the ecoregions map", path);
            //    throw new System.ApplicationException(mesg);
            //}

            using (map)
            {
                DoublePixel pixel = map.BufferPixel;
                foreach (Site site in PlugIn.ModelCore.Landscape.AllSites)
                {
                    map.ReadBufferPixel();
                    double mapValue = pixel.MapCode.Value;
                    if (site.IsActive)
                    {
                        if (mapValue < 0.0001 || mapValue > 0.75)
                            throw new InputValueException(mapValue.ToString(),
                                                          "Soil field capacity value {0} is not between {1:0.0} and {2:0.00}. Site_Row={3:0}, Site_Column={4:0}",
                                                          mapValue, 0.00001, 0.75, site.Location.Row, site.Location.Column);
                        SiteVars.SoilFieldCapacity[site] = mapValue;
                    }
                }
            }
        }
        //---------------------------------------------------------------------

        public static void ReadWiltingPointMap(string path)
        {
            IInputRaster<DoublePixel> map = MakeDoubleMap(path); 

            //try
            //{
            //    map = PlugIn.ModelCore.OpenRaster<DoublePixel>(path);
            //}
            //catch (FileNotFoundException)
            //{
            //    string mesg = string.Format("Error: The file {0} does not exist", path);
            //    throw new System.ApplicationException(mesg);
            //}

            //if (map.Dimensions != PlugIn.ModelCore.Landscape.Dimensions)
            //{
            //    string mesg = string.Format("Error: The input map {0} does not have the same dimension (row, column) as the ecoregions map", path);
            //    throw new System.ApplicationException(mesg);
            //}

            using (map)
            {
                DoublePixel pixel = map.BufferPixel;
                foreach (Site site in PlugIn.ModelCore.Landscape.AllSites)
                {
                    map.ReadBufferPixel();
                    double mapValue = pixel.MapCode.Value;
                    if (site.IsActive)
                    {
                        if (mapValue < 0.0 || mapValue > 0.75)
                            throw new InputValueException(mapValue.ToString(),
                                                          "Soil wilting point value {0} is not between {1:0.0} and {2:0.00}. Site_Row={3:0}, Site_Column={4:0}",
                                                          mapValue, 0.0, 0.75, site.Location.Row, site.Location.Column);
                        if (mapValue > SiteVars.SoilFieldCapacity[site])
                            throw new InputValueException(mapValue.ToString(),
                                                          "Wilting Point {0} is greater than field capacity {1:0.0}.  Site_Row={3:0}, Site_Column={4:0}",
                                                          mapValue, SiteVars.SoilFieldCapacity[site], site.Location.Row, site.Location.Column);
                        SiteVars.SoilWiltingPoint[site] = mapValue;
                    }
                }
            }
        }
        //---------------------------------------------------------------------

        public static void ReadPercentSandMap(string path)
        {
            IInputRaster<DoublePixel> map = MakeDoubleMap(path); 

            //try
            //{
            //    map = PlugIn.ModelCore.OpenRaster<DoublePixel>(path);
            //}
            //catch (FileNotFoundException)
            //{
            //    string mesg = string.Format("Error: The file {0} does not exist", path);
            //    throw new System.ApplicationException(mesg);
            //}

            //if (map.Dimensions != PlugIn.ModelCore.Landscape.Dimensions)
            //{
            //    string mesg = string.Format("Error: The input map {0} does not have the same dimension (row, column) as the ecoregions map", path);
            //    throw new System.ApplicationException(mesg);
            //}

            using (map)
            {
                DoublePixel pixel = map.BufferPixel;
                foreach (Site site in PlugIn.ModelCore.Landscape.AllSites)
                {
                    map.ReadBufferPixel();
                    double mapValue = pixel.MapCode.Value;
                    if (site.IsActive)
                    {
                        if (mapValue < 0.0 || mapValue > 1.0)
                            throw new InputValueException(mapValue.ToString(),
                                                          "Soil percent sand value {0} is not between {1:0.0} and {2:0.0}.  Site_Row={3:0}, Site_Column={4:0}",
                                                          mapValue, 0.0, 1.0, site.Location.Row, site.Location.Column);
                        SiteVars.SoilPercentSand[site] = mapValue;
                    }
                }
            }
        }
        //---------------------------------------------------------------------

        public static void ReadPercentClayMap(string path)
        {
            IInputRaster<DoublePixel> map = MakeDoubleMap(path);

            using (map)
            {
                DoublePixel pixel = map.BufferPixel;
                foreach (Site site in PlugIn.ModelCore.Landscape.AllSites)
                {
                    map.ReadBufferPixel();
                    double mapValue = pixel.MapCode.Value;
                    if (site.IsActive)
                    {
                        if (mapValue < 0.0 || mapValue > 1.0)
                            throw new InputValueException(mapValue.ToString(),
                                                          "Soil percent clay value {0} is not between {1:0.0} and {2:0.0}. Site_Row={3:0}, Site_Column={4:0}",
                                                          mapValue, 0.0, 1.0, site.Location.Row, site.Location.Column);
                        SiteVars.SoilPercentClay[site] = mapValue;
                    }
                }
            }
        }
        //---------------------------------------------------------------------

        public static void ReadSoilBulkDensityMap(string path)
        {
            IInputRaster<DoublePixel> map = MakeDoubleMap(path);

            using (map)
            {
                DoublePixel pixel = map.BufferPixel;
                foreach (Site site in PlugIn.ModelCore.Landscape.AllSites)
                {
                    map.ReadBufferPixel();
                    double mapValue = pixel.MapCode.Value;
                    if (site.IsActive)
                    {
                        if (mapValue < 0.0 || mapValue > 10.0)
                            throw new InputValueException(mapValue.ToString(),
                                                          "Bulk density value {0} is not between {1:0.0} and {2:0.0}. Site_Row={3:0}, Site_Column={4:0}",
                                                          mapValue, 0.0, 10.0, site.Location.Row, site.Location.Column);
                        SiteVars.SoilBulkDensity[site] = mapValue;
                    }
                }
            }
        }

        //---------------------------------------------------------------------

        public static void ReadSoilParticleDensityMap(string path)
        {
            IInputRaster<DoublePixel> map = MakeDoubleMap(path);

            using (map)
            {
                DoublePixel pixel = map.BufferPixel;
                foreach (Site site in PlugIn.ModelCore.Landscape.AllSites)
                {
                    map.ReadBufferPixel();
                    double mapValue = pixel.MapCode.Value;
                    if (site.IsActive)
                    {
                        if (mapValue < 0.0 || mapValue > 10.0)
                            throw new InputValueException(mapValue.ToString(),
                                                          "Particle density value {0} is not between {1:0.0} and {2:0.0}. Site_Row={3:0}, Site_Column={4:0}",
                                                          mapValue, 0.0, 10.0, site.Location.Row, site.Location.Column);
                        SiteVars.SoilParticleDensity[site] = mapValue;
                    }
                }
            }
        }

        public static void ReadSoilCNMaps(string pathSOC, string pathSON)
        {
            IInputRaster<DoublePixel> map = MakeDoubleMap(pathSOC);

            
            //map = MakeDoubleMap(pathSOC);

            using (map)
            {
                DoublePixel pixel = map.BufferPixel;
                foreach (Site site in PlugIn.ModelCore.Landscape.AllSites)
                {
                    map.ReadBufferPixel();
                    double mapValue = pixel.MapCode.Value;
                    if (site.IsActive)
                    {
                        if (mapValue <= 1.0 || mapValue > 10000.0)
                            throw new InputValueException(mapValue.ToString(),
                                                          "SOC value {0} is not between {1:0.0} and {2:0.0}. Site_Row={3:0}, Site_Column={4:0}",
                                                          mapValue, 1.0, 10000.0, site.Location.Row, site.Location.Column);
                        SiteVars.OHorizon[site].Carbon = mapValue;
                    }
                }
            }

            map = MakeDoubleMap(pathSON);

            using (map)
            {
                DoublePixel pixel = map.BufferPixel;
                foreach (Site site in PlugIn.ModelCore.Landscape.AllSites)
                {
                    map.ReadBufferPixel();
                    double mapValue = pixel.MapCode.Value;
                    if (site.IsActive)
                    {
                        if (mapValue <= 0.0 || mapValue > 500.0)
                            throw new InputValueException(mapValue.ToString(),
                                                          "SON value {0} is not between {1:0.0} and {2:0.0}. Site_Row={3:0}, Site_Column={4:0}",
                                                          mapValue, 0.0, 500.0, site.Location.Row, site.Location.Column);
                        SiteVars.OHorizon[site].Nitrogen = mapValue;
                    }
                }
            }

            // Initialize other components:
            foreach (Site site in PlugIn.ModelCore.Landscape.AllSites)
            {
                SiteVars.OHorizon[site].EnzymaticConcentration = 0.001;
                SiteVars.OHorizon[site].MicrobialCarbon = 0.00001;
                SiteVars.OHorizon[site].MicrobialNitrogen = 0.000001;
            }
        }
        ////---------------------------------------------------------------------
        private static IInputRaster<DoublePixel> MakeDoubleMap(string path)
        {
            PlugIn.ModelCore.UI.WriteLine("  Read in data from {0}", path);

            IInputRaster<DoublePixel> map;

            try
            {
                map = PlugIn.ModelCore.OpenRaster<DoublePixel>(path);
            }
            catch (FileNotFoundException)
            {
                string mesg = string.Format("Error: The file {0} does not exist", path);
                throw new System.ApplicationException(mesg);
            }

            if (map.Dimensions != PlugIn.ModelCore.Landscape.Dimensions)
            {
                string mesg = string.Format("Error: The input map {0} does not have the same dimension (row, column) as the scenario ecoregions map", path);
                throw new System.ApplicationException(mesg);
            }

            return map;
        }
        //---------------------------------------------------------------------

        public static void ReadDeadWoodMaps(string surfacePath, string soilPath)
        {
            IInputRaster<DoublePixel> map = MakeDoubleMap(surfacePath);

            using (map)
            {
                DoublePixel pixel = map.BufferPixel;
                foreach (Site site in PlugIn.ModelCore.Landscape.AllSites)
                {
                    map.ReadBufferPixel();
                    double mapValue = pixel.MapCode.Value;
                    if (site.IsActive)
                    {
                        if (mapValue < 0.0 || mapValue > 50000.0)
                            throw new InputValueException(mapValue.ToString(),
                                                          "SurfDeadWood value {0} is not between {1:0.0} and {2:0.0}. Site_Row={3:0}, Site_Column={4:0}",
                                                          mapValue, 0.0, 50000.0, site.Location.Row, site.Location.Column);
                        SiteVars.SurfaceDeadWood[site].Carbon = mapValue * 0.47;
                        SiteVars.SurfaceDeadWood[site].Nitrogen = mapValue * 0.47 / 200.0;  // 200 is a generic wood CN ratio
                        SiteVars.SoilStructural[site].Carbon = SiteVars.SurfaceDeadWood[site].Carbon * 0.85 * PlugIn.Parameters.InitialFineFuels;
                        SiteVars.SoilStructural[site].Nitrogen = SiteVars.SoilStructural[site].Carbon / OtherData.StructuralCN;
                        SiteVars.SoilMetabolic[site].Carbon = SiteVars.SurfaceDeadWood[site].Carbon * 0.15 * PlugIn.Parameters.InitialFineFuels;
                        SiteVars.SoilMetabolic[site].Nitrogen = SiteVars.SoilMetabolic[site].Carbon / 10.0;  // a generic metabolic CN ratio

                    }
                }
            }

            map = MakeDoubleMap(soilPath);

            // Soil Dead Wood = Dead Roots

            using (map)
            {
                DoublePixel pixel = map.BufferPixel;
                foreach (Site site in PlugIn.ModelCore.Landscape.AllSites)
                {
                    map.ReadBufferPixel();
                    double mapValue = pixel.MapCode.Value;
                    if (site.IsActive)
                    {
                        if (mapValue < 0.0 || mapValue > 50000.0)
                            throw new InputValueException(mapValue.ToString(),
                                                          "SoilDeadWood value {0} is not between {1:0.0} and {2:0.0}. Site_Row={3:0}, Site_Column={4:0}",
                                                          mapValue, 0.0, 50000.0, site.Location.Row, site.Location.Column);
                        SiteVars.SoilDeadWood[site].Carbon = mapValue * 0.47;
                        SiteVars.SoilDeadWood[site].Nitrogen = mapValue * 0.47 / 200.0;  // 200 is a generic wood CN ratio
                    }
                }
            }
        }
        //---------------------------------------------------------------------
        //---------------------------------------------------------------------

        //public static void ReadDoubleMap(string path, ISiteVar<double> siteVar)
        //{
        //    IInputRaster<DoublePixel> map;

        //    try
        //    {
        //        map = PlugIn.ModelCore.OpenRaster<DoublePixel>(path);
        //    }
        //    catch (FileNotFoundException)
        //    {
        //        string messege = string.Format("Error: The file {0} does not exist", path);
        //        throw new System.ApplicationException(messege);
        //    }

        //    if (map.Dimensions != PlugIn.ModelCore.Landscape.Dimensions)
        //    {
        //        string messege = string.Format("Error: The input map {0} does not have the same dimension (row, column) as the ecoregions map", path);
        //        throw new System.ApplicationException(messege);
        //    }

        //    using (map)
        //    {
        //        DoublePixel pixel = map.BufferPixel;
        //        foreach (Site site in PlugIn.ModelCore.Landscape.AllSites)
        //        {
        //            map.ReadBufferPixel();
        //            double mapCode = pixel.MapCode.Value;

        //            if (site.IsActive)
        //            {
        //                siteVar[site] = mapCode;
        //            }
                //}
            //}
        //}
    }
}
