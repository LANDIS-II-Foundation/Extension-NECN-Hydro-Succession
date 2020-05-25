# Extension-NECN-Succession

The next generation NECN, v7.x

Changes to this extension are governed by the [**Repository Rules**](https://sites.google.com/site/landismodel/developers) from the Technical Advisory Committee.

Robert Scheller is the Maintainer of this repository.

# Primary Differences with NECN v6

NECN v7 will feauture simplified soil respiration algorithms.  A new OHorizon uses the respriation function outlined by the DAMM-MnCIP model (Abramoff and others insert REF).  

What was SOM2 has been recoded as the MineralSoil layer.  This layer follows the general principle of Century respiration.

The surface litter and soil litter (both containing lignin and metabolic components) and wood pools have also been retained and follow the Century decay functions.

# Specific Differences with NECN v6

Note: Listed below are the substantive changes.  Other minor changes (names of variables, debug logging statements) located throughout the code.

DAMMLayer.cs: Contains the OHorizon properties and respiration function. Should be renamed OrganicHorizonLayer.cs L362: Following Century, a fraction goes to MineralSoil (a user input).

Layer.cs:  DecomposeLignin:  Lignin C goes first to MineralSoil and the remainder to OHorizon.  This follows the Century order-of-operations.  This order could be reversed by reconsidering the meaning of OtherData.LigninRespirationRate.

Layer.cs:  DecomposeMetabolic:  ML added flow to DOC, L273-280.  RMS:  Divide between OHorizon and MineralSoil depending on soil depth and assumed OHorizon depth (10 cm), L298.  Allocated to respective pools L304-309.

Main.cs:  Added OHorizon decompose, L108.  Order goes from top -> down.

MineralSoilLayer.cs (nee SoilLayer.cs):  The Decompose function has been radically simplified reflecting the elimination of the SOM1 and SOM3 layers.  L44:  Assume decomposition always possible.  L49: Following Century, define a fraction that goes to CO2 (this is now a user input rather than a fixed constant).  UNRESOLVED:  Where does the remainder go (L50)?

Outputs.cs:  L66: Primary Log File restructured to reflect the changes above.

PlugIn.cs:  L101:  Also initialize MineralSoil layer.




