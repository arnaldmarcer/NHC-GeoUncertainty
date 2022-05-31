## Scripts for:  

**Uncertainty matters: ascertaining where specimens in natural history collections come from.**  
Marcer A. (1,2), Chapman A.D. (3), Wieczorek J.R. (4), Picó F.X. (5), Uribe F. (6), Waller J. (7) & Ariño A.H. (8,9)

(1) CREAF, E08193 Bellaterra (Cerdanyola del Vallès), Catalonia, Spain  
(2) Universitat Autònoma de Barcelona, E08193 Bellaterra (Cerdanyola del Vallès), Catalonia, Spain  
(3) Australian Biodiversity Information Services, Melbourne, Victoria, Australia  
(4) University of California, Berkeley, California, USA  
(5) Estación Biológica de Doñana (EBD), Consejo Superior de Investigaciones Científicas (CSIC), Sevilla, Spain  
(6) Museu de Ciències Naturals, Barcelona, Catalonia, Spain  
(7) GBIF, Copenhagen, Denmark  
(8) Institute for Biodiversity and Environmental Research (BIOMA) and DATAI, Universidad de Navarra, Pamplona, Spain  
(9) Museo de Ciencias de la Universidad de Navarra, Pamplona, Spain  

*Ecography, 2022*

---
### Order of execution of the scripts  
1. scripts/bash/gbif_data.sh  
2. scripts/R/non-geo-dataset-preparation.R  
3. scripts/R/non-geo-results-preparation.R  
4. scripts/R/geo-dataset-preparation.R  
5. scripts/R/geo-results-preparation.R  
6. scripts/R/tree_cover.R  
7. scripts/R/sdm.R  
8. scripts/R/suitability-scores.R  
9. scripts/R/figures-and-tables.R  
---

### Main folder structure  

- NHC-GeoUncertainty  
    - data  
        - raw 
            - gbif  
            - global_forest_cover  
            - worldlim  
        - processed  
            - gbif  
            - predictors  
                - australia  
                - mexico_to_argentina  
                - northern_canada_and_greenland  
    - manuscript  
        - figures  
        - tables  
    - outputs  
        - maxent  
            - eucalyptus_gongylocarpa-uncertainty  
            - guazuma_ulmifolia-uncertainty  
            - rhododendron_groenlandicum-uncertainty  
    - scripts  
        - bash  
        - R  
    - tmp  
---
