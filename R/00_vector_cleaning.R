# Cleaning and standardisation of the land use/land cover vector (Vercelli sector).
# The raw vector (3307 polygons, 53 classes) carries up to 20 DESCRIZION variants
# per CLASSE code (typos, question marks, notes). This script maps each code to one
# clean description, checks that no polygon is lost or duplicated, and saves a
# _CLEAN.gpkg (the original file is left untouched). Geometry and the CLASSE field
# are never modified.
# Note: class 3111 merges two formations that share the same code (Salix alba
# willows and mesohygrophilous riparian stands); to keep them distinct, add a
# dedicated field in QGIS before running this script.

# Packages ----
library(terra)   # read/write geographic vectors (.gpkg)
library(dplyr)   # manipulate the attribute table as a data frame

# Folders ----
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "
     Set the path at the top of the script, or the POPARK_DATA variable.")
dir_output <- file.path(DIR_BASE, "output_spettrale")

file_input  <- file.path(DIR_BASE, "usosuolo_2019_layout_def.shp")        # raw input (never modified)
file_output <- file.path(DIR_BASE, "usosuolo_2019_layout_def_CLEAN.gpkg") # clean output
file_verif  <- file.path(dir_output, "verifica_pulizia.csv")             # verification CSV

# Load the original vector ----
cat("Loading original vector...\n")
usosuolo <- vect(file_input)

# store the original polygon count to verify none is lost or duplicated later
n_originale <- nrow(usosuolo)
cat("Total polygons:", n_originale, "\n")   # expected: 3307

# work on the attribute table only; the geometry is never touched
df <- as.data.frame(usosuolo)

# stop with a clear message if the required fields are missing
if (!"CLASSE" %in% names(df)) {
  stop("ERROR: field 'CLASSE' not found in the layer. ",
       "Available fields: ", paste(names(df), collapse = ", "))
}
if (!"DESCRIZION" %in% names(df)) {
  stop("ERROR: field 'DESCRIZION' not found in the layer. ",
       "Available fields: ", paste(names(df), collapse = ", "))
}

# before cleaning we expect 53 classes and 161 CLASSE+DESCRIZION variants
cat("Unique classes (CLASSE) before cleaning:", length(unique(df$CLASSE)), "\n")
cat("CLASSE+DESCRIZION variants before cleaning:",
    nrow(distinct(df, CLASSE, DESCRIZION)), "\n\n")

# Standardisation table (lookup) ----
# Maps each CLASSE code to a single canonical description (most frequent, spelling
# corrected, no question marks or notes). Exactly 53 rows, one per CLASSE code.
# Do not change the CLASSE column: these are the original codes and must match.
lookup <- data.frame(
  CLASSE = c(
    # artificial areas (100-199)
     111,  112,  121,  122,  142,
    # agricultural areas (200-299)
     213,  221,  222,  224,  225,  226,
    # open spaces with sparse vegetation (300-399, first group)
     332,  335,
    # quarries, dumps, degraded areas (1300-1499)
    1311, 1312, 1313, 1314, 1315,
    1321, 1322, 1421, 1422,
    # arable land and crops (2100-2399)
    2111, 2112, 2113, 2121, 2311,
    # woodlands and riparian formations (3100-3199)
    3111, 3112, 3113, 3114, 3115, 3116, 3117,
    # herbaceous and shrub vegetation (3200-3299)
    3212,
    # shrub and transitional formations (3240-3299)
    3241, 3242, 3243, 3244, 3246,
    # sands and gravel bars (3300-3399)
    3311, 3312,
    # wetlands (4100-4299)
    4113, 4241, 4242, 4243, 4244,
    # water bodies (5100-5199)
    5111, 5112, 5113, 5121, 5122, 5123
  ),
  DESC_CLEAN = c(
    # artificial areas
    "Tessuto urbano continuo",           # 111: 16 polygons, 1 variant - ok
    "Tessuto urbano discontinuo",        # 112: 114 polygons, 1 variant - ok
    "Aree industriali e commerciali",    # 121: 11 polygons, 2 variants ('e' vs 'o')
    "Reti stradali e ferroviarie e spazi accessori", # 122: 14 polygons, 2 variants (typo 'ferroviari')
    "Aree sportive e ricreative",        # 142: 1 polygon, 1 variant - ok

    # agricultural areas
    "Risaie",                            # 213: 31 polygons, 1 variant - ok
    "Vigneti",                           # 221: 17 polygons, 1 variant - ok
    "Frutteti e frutti minori",          # 222: 82 polygons, 2 variants (one with 'noccioleti??')
    "Pioppeti",                          # 224: 322 polygons, 5 variants (typo and ?)
    "Altri impianti di arboricoltura da legno", # 225: 79 polygons, 3 variants (typo)
    "Terreni non in coltivazione",       # 226: 198 polygons, 5 variants (typo and ?)

    # open spaces
    "Rocce nude",                        # 332: 14 polygons, 2 variants ('(Calanchi)' removed)
    "Aree interessate da processi erosionali e depositi diffusi", # 335: 42 polygons, 12 variants!

    # quarries and extraction areas
    "Piazzali di cava",                  # 1311: 9 polygons, 1 variant - ok
    "Specchi d'acqua di cave attive",    # 1312: 1 polygon, 1 variant - ok
    "Specchi d'acqua di cave non rinaturalizzate",  # 1313: 3 polygons, 3 variants (typo)
    "Spazi seminaturali all'interno di cave attive", # 1314: 1 polygon, 1 variant - ok
    "Aree interessate da scotico superficiale",      # 1315: 7 polygons, 3 variants (typo)

    # dumps and degraded areas
    "Discariche",                        # 1321: 3 polygons, 1 variant - ok
    "Aree degradate",                    # 1322: 9 polygons, 2 variants ('Area' vs 'Aree')
    "Aree con baracche",                 # 1421: 33 polygons, 2 variants ('Area' vs 'Aree')
    "Altre aree ricreative e sportive",  # 1422: 7 polygons, 2 variants (word order)

    # arable land and crops
    "Orticoltura in pieno campo",        # 2111: 15 polygons, 1 variant - ok
    "Seminativi in aree non irrigue",    # 2112: 264 polygons, 2 variants (one with 'Arundo donax???')
    "Serre",                             # 2113: 9 polygons, 1 variant - ok
    "Seminativi in aree irrigue",        # 2121: 26 polygons, 2 variants (double space)
    "Prati stabili",                     # 2311: 32 polygons, 4 variants (typo and ?)

    # woodlands and riparian formations
    # 3111 groups 'Saliceti arborei a Salix alba' (215) and 'Formazioni arboree
    # riparie mesoigrofile' (111) under the same code
    "Saliceti arborei e formazioni arboree riparie", # 3111: 336 polygons, 7 variants
    "Alneti e formazioni arboree igrofile delle lanche", # 3112: 38 polygons, 3 variants (?)
    "Formazioni arboree planiziali",     # 3113: 33 polygons, 2 variants (???)
    "Boschi misti collinari",            # 3114: 33 polygons, 3 variants (notes and ?)
    "Robinieti",                         # 3115: 130 polygons, 7 variants (from 1 to 4 ?)
    "Filari arborei",                    # 3116: 38 polygons, 2 variants (??)
    "Imboschimenti a conifere",          # 3117: 1 polygon, 1 variant - ok

    # tall herbaceous vegetation (3212 has 20 variants, the worst fragmentation)
    "Popolamenti alto-erbacei",          # 3212: 272 polygons, 20 variants (typos, capitals, spaces)

    # shrub and transitional formations
    "Formazioni erbaceo-arbustive xerofile stabili", # 3241: 43 polygons, 7 variants (?)
    "Saliceti arbustivi",                # 3242: 84 polygons, 7 variants (?, Arundo, reforestation)
    "Formazioni arbustive di ricolonizzazione", # 3243: 284 polygons, 5 variants (typo)
    "Imboschimenti",                     # 3244: 54 polygons, 2 variants (capital)
    "Siepi",                             # 3246: 56 polygons, 2 variants ('Siepi/fasce rispetto')

    # sands, gravel, riverbed
    "Sabbie e ghiaioni",                 # 3311: 123 polygons, 1 variant - ok
    "Greto",                             # 3312: 131 polygons, 3 variants ('Greti', 'con Arundo donax')

    # wetlands
    "Formazioni erbaceo-arbustive a dominanza di igrofite (lanche)", # 4113: 1 polygon
    "Lanche",                            # 4241: 62 polygons, 1 variant - ok
    "Formazioni erbaceo-arbustive a dominanza di igrofite", # 4242: 45 polygons, 7 variants
    "Popolamenti vegetali acquatici",    # 4243: 44 polygons, 1 variant - ok
    "Formazioni a elofite",              # 4244: 4 polygons, 1 variant - ok

    # water bodies
    "Alveo principale",                  # 5111: 15 polygons, 1664 ha - the main river
    "Bracci secondari",                  # 5112: 23 polygons, 1 variant - ok
    "Canali",                            # 5113: 44 polygons, 4 variants (with/without 'pertinenze')
    "Canali e solchi erosionali",        # 5121: 12 polygons, 4 variants (typo 'eorsionali')
    "Specchi d'acqua da cave rinaturalizzate",  # 5122: 13 polygons, 3 variants (typo and ?)
    "Specchi d'acqua artificiali"        # 5123: 13 polygons, 1 variant - ok
  ),
  stringsAsFactors = FALSE   # keep strings as strings (no factor conversion)
)

# Check that the lookup is complete ----
# if the vector has a code missing from the lookup, stop before any modification
classi_vett     <- sort(unique(df$CLASSE))
classi_lookup   <- sort(lookup$CLASSE)
classi_mancanti <- setdiff(classi_vett, classi_lookup)

if (length(classi_mancanti) > 0) {
  cat("! WARNING: these CLASSE codes are in the vector but NOT in the lookup:\n")
  print(classi_mancanti)
  stop("Script stopped: incomplete lookup. Add these codes before proceeding.")
} else {
  cat("OK all", length(classi_vett), "CLASSE codes of the vector are in the lookup.\n\n")
}

# Apply the standardisation ----
# drop the old DESCRIZION (with typos), join the clean description from the lookup,
# then rename it to DESCRIZION so the attribute table structure stays identical
cat("Applying description standardisation...\n")
df_clean <- df %>%
  select(-DESCRIZION) %>%
  left_join(lookup, by = "CLASSE") %>%
  rename(DESCRIZION = DESC_CLEAN)

# check 1: polygon count must be unchanged (a change means duplicated lookup codes)
if (nrow(df_clean) != n_originale) {
  stop("CRITICAL ERROR: the polygon count changed during cleaning!\n",
       "Before: ", n_originale, " | After: ", nrow(df_clean), "\n",
       "Check the lookup for duplicated CLASSE codes.")
}
cat("OK polygon count unchanged:", nrow(df_clean), "\n")

# check 2: after cleaning, unique CLASSE+DESCRIZION combinations must equal classes
n_varianti_dopo <- nrow(distinct(df_clean, CLASSE, DESCRIZION))
n_classi_dopo   <- length(unique(df_clean$CLASSE))
cat("Unique classes after cleaning:", n_classi_dopo, "\n")
cat("CLASSE+DESCRIZION variants after cleaning:", n_varianti_dopo, "\n")
if (n_varianti_dopo != n_classi_dopo) {
  stop("ERROR: some classes still have more than one description. Check the lookup.")
} else {
  cat("OK every CLASSE code has exactly one canonical description.\n\n")
}

# Rebuild the vector and save ----
# replace the attribute table with the clean one; geometry is not touched
cat("Rebuilding the vector with the clean attribute table...\n")
values(usosuolo) <- df_clean

cat("\nFirst 5 rows of the clean vector (visual check):\n")
print(head(as.data.frame(usosuolo)[, c("CLASSE", "DESCRIZION", "AREA_HA")], 5))

# never overwrites the original (different name)
cat("\nSaving the clean vector...\n")
writeVector(usosuolo, file_output, overwrite = TRUE)
cat("OK clean vector saved to:\n ", file_output, "\n")

# Final summary and verification CSV ----
# one row per class with polygon count and total area (hectares)
cat("\nWriting the final summary CSV...\n")
riepilogo_finale <- as.data.frame(usosuolo) %>%
  group_by(CLASSE, DESCRIZION) %>%
  summarise(
    n_poligoni = n(),
    area_ha    = round(sum(AREA_HA, na.rm = TRUE), 2),
    .groups    = "drop"
  ) %>%
  arrange(CLASSE)

# ';' separator (European format, compatible with Italian Excel)
write.csv2(riepilogo_finale, file_verif, row.names = FALSE)

cat("\nClass summary after cleaning:\n")
cat(sprintf("%-8s  %-6s  %10s  %s\n", "CLASSE", "N_pol", "Area_ha", "Description"))
for (i in seq_len(nrow(riepilogo_finale))) {
  r <- riepilogo_finale[i, ]
  cat(sprintf("%-8s  %-6s  %10.1f  %s\n",
              r$CLASSE, r$n_poligoni, r$area_ha, substr(r$DESCRIZION, 1, 48)))
}
cat("Total classes:  ", nrow(riepilogo_finale), "\n")
cat("Total polygons: ", sum(riepilogo_finale$n_poligoni), "\n")   # expected: 3307
cat("Total area (ha):", round(sum(riepilogo_finale$area_ha), 1), "\n")
cat("OK verification CSV saved to:\n ", file_verif, "\n")

# Final check in QGIS ----
# open the OUTPUT file usosuolo_2019_layout_def_CLEAN.gpkg (the original .shp is
# not modified) and spot-check the attribute table:
#   CLASSE = 335  -> all polygons share the same DESCRIZION
#   CLASSE = 224  -> all read "Pioppeti"
#   CLASSE = 3212 -> all read "Popolamenti alto-erbacei"
