# GEDSpy
#### GEDSpy is the python library for gene list enrichment with genes ontology, pathways and potential drugs

<p align="right">
<img  src="https://github.com/jkubis96/GEDSpy/blob/main/fig/logo_jbs.PNG?raw=true" alt="drawing" width="250" />
</p>


### Author: Jakub Kubiś 

<div align="left">
 Institute of Bioorganic Chemistry<br />
 Polish Academy of Sciences<br />
 Department of Molecular Neurobiology<br />
</div>


<br />
<br />

<div align="justify"> </div>

<br />

<p align="center">
<img  src="https://github.com/jkubis96/GEDSpy/blob/main/fig/pathways_pathway.png?raw=true" alt="drawing" width="600" />
</p>

##### Figure 1 Single-cell sequencing in DropSeq technology A) libraries preparing  B) sequencing and analysis  [Created in BioRender]

<br />


<div align="justify"> </div>


```html

<html>
<head>
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/vis/4.16.1/vis.css" type="text/css" />
<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/vis/4.16.1/vis-network.min.js"> </script>
<center>
<h1></h1>
</center>

<!-- <link rel="stylesheet" href="../node_modules/vis/dist/vis.min.css" type="text/css" />
<script type="text/javascript" src="../node_modules/vis/dist/vis.js"> </script>-->

<style type="text/css">

        #mynetwork {
            width: 1000px;
            height: 800px;
            background-color: #ffffff;
            border: 1px solid lightgray;
            position: relative;
            float: left;
        }

        

        
        #config {
            float: left;
            width: 400px;
            height: 600px;
        }
        

        
</style>

</head>

<body>
<div id = "mynetwork"></div>


<div id = "config"></div>

<script type="text/javascript">

    // initialize global variables.
    var edges;
    var nodes;
    var network; 
    var container;
    var options, data;

    
    // This method is responsible for drawing the graph, returns the drawn network
    function drawGraph() {
        var container = document.getElementById('mynetwork');
        
        

        // parsing and collecting nodes and edges from the python
        nodes = new vis.DataSet([{"color": "red", "id": "EGFR", "label": "EGFR", "shape": "dot", "size": 37}, {"color": "red", "id": "CHEK2", "label": "CHEK2", "shape": "dot", "size": 30}, {"color": "red", "id": "DUSP9", "label": "DUSP9", "shape": "dot", "size": 13}, {"color": "red", "id": "ESR1", "label": "ESR1", "shape": "dot", "size": 34}, {"color": "red", "id": "FGFR2", "label": "FGFR2", "shape": "dot", "size": 25}, {"color": "red", "id": "GNAI2", "label": "GNAI2", "shape": "dot", "size": 29}, {"color": "red", "id": "GNB2", "label": "GNB2", "shape": "dot", "size": 27}, {"color": "red", "id": "GNB3", "label": "GNB3", "shape": "dot", "size": 26}, {"color": "red", "id": "GPNMB", "label": "GPNMB", "shape": "dot", "size": 6}, {"color": "red", "id": "GRB10", "label": "GRB10", "shape": "dot", "size": 13}, {"color": "red", "id": "HSP90AB1", "label": "HSP90AB1", "shape": "dot", "size": 34}, {"color": "red", "id": "IGF2R", "label": "IGF2R", "shape": "dot", "size": 19}, {"color": "red", "id": "IL1R1", "label": "IL1R1", "shape": "dot", "size": 6}, {"color": "red", "id": "IRAK1", "label": "IRAK1", "shape": "dot", "size": 34}, {"color": "red", "id": "ITGA1", "label": "ITGA1", "shape": "dot", "size": 16}, {"color": "red", "id": "ITPR1", "label": "ITPR1", "shape": "dot", "size": 26}, {"color": "red", "id": "ITSN1", "label": "ITSN1", "shape": "dot", "size": 17}, {"color": "red", "id": "JAK1", "label": "JAK1", "shape": "dot", "size": 28}, {"color": "red", "id": "KALRN", "label": "KALRN", "shape": "dot", "size": 21}, {"color": "red", "id": "KRAS", "label": "KRAS", "shape": "dot", "size": 28}, {"color": "red", "id": "KSR1", "label": "KSR1", "shape": "dot", "size": 29}, {"color": "red", "id": "LAMA5", "label": "LAMA5", "shape": "dot", "size": 16}, {"color": "red", "id": "LAMB2", "label": "LAMB2", "shape": "dot", "size": 17}, {"color": "red", "id": "MAP2K7", "label": "MAP2K7", "shape": "dot", "size": 29}, {"color": "red", "id": "MAPK13", "label": "MAPK13", "shape": "dot", "size": 32}, {"color": "red", "id": "MUSK", "label": "MUSK", "shape": "dot", "size": 6}, {"color": "red", "id": "MYO9B", "label": "MYO9B", "shape": "dot", "size": 16}, {"color": "red", "id": "NCSTN", "label": "NCSTN", "shape": "dot", "size": 23}, {"color": "red", "id": "NUP210", "label": "NUP210", "shape": "dot", "size": 6}, {"color": "red", "id": "PIM1", "label": "PIM1", "shape": "dot", "size": 28}, {"color": "red", "id": "PPARG", "label": "PPARG", "shape": "dot", "size": 32}, {"color": "red", "id": "PPP2R1A", "label": "PPP2R1A", "shape": "dot", "size": 30}, {"color": "red", "id": "PRKCE", "label": "PRKCE", "shape": "dot", "size": 35}, {"color": "red", "id": "PTK2B", "label": "PTK2B", "shape": "dot", "size": 31}, {"color": "red", "id": "PTPRK", "label": "PTPRK", "shape": "dot", "size": 6}, {"color": "red", "id": "RARA", "label": "RARA", "shape": "dot", "size": 29}, {"color": "red", "id": "ROCK2", "label": "ROCK2", "shape": "dot", "size": 32}, {"color": "red", "id": "SDC1", "label": "SDC1", "shape": "dot", "size": 17}, {"color": "red", "id": "SLC8A3", "label": "SLC8A3", "shape": "dot", "size": 6}, {"color": "red", "id": "SMARCA2", "label": "SMARCA2", "shape": "dot", "size": 21}, {"color": "red", "id": "STAT1", "label": "STAT1", "shape": "dot", "size": 33}, {"color": "red", "id": "TNFSF14", "label": "TNFSF14", "shape": "dot", "size": 6}, {"color": "red", "id": "TNK2", "label": "TNK2", "shape": "dot", "size": 27}, {"color": "red", "id": "TRIB3", "label": "TRIB3", "shape": "dot", "size": 27}, {"color": "red", "id": "CAMK1G", "label": "CAMK1G", "shape": "dot", "size": 10}, {"color": "red", "id": "CDK18", "label": "CDK18", "shape": "dot", "size": 16}, {"color": "red", "id": "CDK19", "label": "CDK19", "shape": "dot", "size": 13}, {"color": "red", "id": "MAP4K2", "label": "MAP4K2", "shape": "dot", "size": 17}, {"color": "red", "id": "MAP4K3", "label": "MAP4K3", "shape": "dot", "size": 16}, {"color": "red", "id": "MAST2", "label": "MAST2", "shape": "dot", "size": 17}, {"color": "red", "id": "SNRK", "label": "SNRK", "shape": "dot", "size": 17}, {"color": "red", "id": "TLK2", "label": "TLK2", "shape": "dot", "size": 17}, {"color": "red", "id": "DBNL", "label": "DBNL", "shape": "dot", "size": 6}, {"color": "red", "id": "PFKL", "label": "PFKL", "shape": "dot", "size": 6}, {"color": "red", "id": "DUSP10", "label": "DUSP10", "shape": "dot", "size": 10}, {"color": "red", "id": "PTPN23", "label": "PTPN23", "shape": "dot", "size": 6}, {"color": "red", "id": "EGR2", "label": "EGR2", "shape": "dot", "size": 21}, {"color": "red", "id": "KDM3A", "label": "KDM3A", "shape": "dot", "size": 19}, {"color": "red", "id": "NCOA6", "label": "NCOA6", "shape": "dot", "size": 16}, {"color": "red", "id": "PEX14", "label": "PEX14", "shape": "dot", "size": 6}, {"color": "red", "id": "PPARGC1B", "label": "PPARGC1B", "shape": "dot", "size": 17}, {"color": "red", "id": "RNF14", "label": "RNF14", "shape": "dot", "size": 13}, {"color": "red", "id": "RREB1", "label": "RREB1", "shape": "dot", "size": 13}, {"color": "red", "id": "SRF", "label": "SRF", "shape": "dot", "size": 19}, {"color": "red", "id": "SPR", "label": "SPR", "shape": "dot", "size": 6}, {"color": "red", "id": "CD46", "label": "CD46", "shape": "dot", "size": 6}, {"color": "red", "id": "SERPINA1", "label": "SERPINA1", "shape": "dot", "size": 6}, {"color": "red", "id": "COL6A3", "label": "COL6A3", "shape": "dot", "size": 10}, {"color": "red", "id": "ITGA7", "label": "ITGA7", "shape": "dot", "size": 6}, {"color": "red", "id": "SLCO1A2", "label": "SLCO1A2", "shape": "dot", "size": 6}, {"color": "red", "id": "SLCO1B3", "label": "SLCO1B3", "shape": "dot", "size": 6}]);
        edges = new vis.DataSet([{"from": "EGFR", "to": "CHEK2", "weight": 24.849066497880003}, {"from": "EGFR", "to": "DUSP9", "weight": 21.972245773362197}, {"from": "EGFR", "to": "ESR1", "weight": 31.780538303479457}, {"from": "EGFR", "to": "FGFR2", "weight": 30.91042453358316}, {"from": "EGFR", "to": "GNAI2", "weight": 27.725887222397812}, {"from": "EGFR", "to": "GNB2", "weight": 27.0805020110221}, {"from": "EGFR", "to": "GNB3", "weight": 23.978952727983707}, {"from": "EGFR", "to": "GPNMB", "weight": 21.972245773362197}, {"from": "EGFR", "to": "GRB10", "weight": 23.02585092994046}, {"from": "EGFR", "to": "HSP90AB1", "weight": 30.91042453358316}, {"from": "EGFR", "to": "IGF2R", "weight": 27.0805020110221}, {"from": "EGFR", "to": "IL1R1", "weight": 23.02585092994046}, {"from": "EGFR", "to": "IRAK1", "weight": 27.725887222397812}, {"from": "EGFR", "to": "ITGA1", "weight": 24.849066497880003}, {"from": "EGFR", "to": "ITPR1", "weight": 27.725887222397812}, {"from": "EGFR", "to": "ITSN1", "weight": 23.978952727983707}, {"from": "EGFR", "to": "JAK1", "weight": 29.444389791664403}, {"from": "EGFR", "to": "KALRN", "weight": 21.972245773362197}, {"from": "EGFR", "to": "KRAS", "weight": 30.44522437723423}, {"from": "EGFR", "to": "KSR1", "weight": 26.390573296152585}, {"from": "EGFR", "to": "LAMA5", "weight": 21.972245773362197}, {"from": "EGFR", "to": "LAMB2", "weight": 23.02585092994046}, {"from": "EGFR", "to": "MAP2K7", "weight": 28.903717578961647}, {"from": "EGFR", "to": "MAPK13", "weight": 24.849066497880003}, {"from": "EGFR", "to": "MUSK", "weight": 23.02585092994046}, {"from": "EGFR", "to": "MYO9B", "weight": 23.02585092994046}, {"from": "EGFR", "to": "NCSTN", "weight": 25.649493574615366}, {"from": "EGFR", "to": "NUP210", "weight": 21.972245773362197}, {"from": "EGFR", "to": "PIM1", "weight": 23.02585092994046}, {"from": "EGFR", "to": "PPARG", "weight": 29.444389791664403}, {"from": "EGFR", "to": "PPP2R1A", "weight": 30.44522437723423}, {"from": "EGFR", "to": "PRKCE", "weight": 27.725887222397812}, {"from": "EGFR", "to": "PTK2B", "weight": 29.957322735539908}, {"from": "EGFR", "to": "PTPRK", "weight": 23.02585092994046}, {"from": "EGFR", "to": "RARA", "weight": 28.903717578961647}, {"from": "EGFR", "to": "ROCK2", "weight": 27.0805020110221}, {"from": "EGFR", "to": "SDC1", "weight": 26.390573296152585}, {"from": "EGFR", "to": "SLC8A3", "weight": 23.02585092994046}, {"from": "EGFR", "to": "SMARCA2", "weight": 21.972245773362197}, {"from": "EGFR", "to": "STAT1", "weight": 30.44522437723423}, {"from": "EGFR", "to": "TNFSF14", "weight": 23.02585092994046}, {"from": "EGFR", "to": "TNK2", "weight": 25.649493574615366}, {"from": "EGFR", "to": "TRIB3", "weight": 23.978952727983707}, {"from": "PRKCE", "to": "CAMK1G", "weight": 21.972245773362197}, {"from": "PRKCE", "to": "CDK18", "weight": 21.972245773362197}, {"from": "PRKCE", "to": "CDK19", "weight": 21.972245773362197}, {"from": "PRKCE", "to": "CHEK2", "weight": 24.849066497880003}, {"from": "PRKCE", "to": "ESR1", "weight": 24.849066497880003}, {"from": "PRKCE", "to": "FGFR2", "weight": 21.972245773362197}, {"from": "PRKCE", "to": "GNAI2", "weight": 26.390573296152585}, {"from": "PRKCE", "to": "GNB2", "weight": 23.978952727983707}, {"from": "PRKCE", "to": "GNB3", "weight": 23.978952727983707}, {"from": "PRKCE", "to": "HSP90AB1", "weight": 23.978952727983707}, {"from": "PRKCE", "to": "IRAK1", "weight": 27.725887222397812}, {"from": "PRKCE", "to": "ITPR1", "weight": 27.725887222397812}, {"from": "PRKCE", "to": "ITSN1", "weight": 23.02585092994046}, {"from": "PRKCE", "to": "JAK1", "weight": 23.02585092994046}, {"from": "PRKCE", "to": "KALRN", "weight": 26.390573296152585}, {"from": "PRKCE", "to": "KSR1", "weight": 27.725887222397812}, {"from": "PRKCE", "to": "MAP2K7", "weight": 26.390573296152585}, {"from": "PRKCE", "to": "MAP4K2", "weight": 23.02585092994046}, {"from": "PRKCE", "to": "MAP4K3", "weight": 21.972245773362197}, {"from": "PRKCE", "to": "MAPK13", "weight": 28.903717578961647}, {"from": "PRKCE", "to": "MAST2", "weight": 23.02585092994046}, {"from": "PRKCE", "to": "MYO9B", "weight": 23.02585092994046}, {"from": "PRKCE", "to": "NCSTN", "weight": 21.972245773362197}, {"from": "PRKCE", "to": "PIM1", "weight": 26.390573296152585}, {"from": "PRKCE", "to": "PPARG", "weight": 24.849066497880003}, {"from": "PRKCE", "to": "PPP2R1A", "weight": 27.0805020110221}, {"from": "PRKCE", "to": "PTK2B", "weight": 26.390573296152585}, {"from": "PRKCE", "to": "RARA", "weight": 21.972245773362197}, {"from": "PRKCE", "to": "ROCK2", "weight": 28.33213344056216}, {"from": "PRKCE", "to": "SNRK", "weight": 23.02585092994046}, {"from": "PRKCE", "to": "STAT1", "weight": 23.978952727983707}, {"from": "PRKCE", "to": "TLK2", "weight": 23.978952727983707}, {"from": "PRKCE", "to": "TNK2", "weight": 25.649493574615366}, {"from": "PRKCE", "to": "TRIB3", "weight": 24.849066497880003}, {"from": "HSP90AB1", "to": "CHEK2", "weight": 21.972245773362197}, {"from": "HSP90AB1", "to": "DBNL", "weight": 21.972245773362197}, {"from": "HSP90AB1", "to": "ESR1", "weight": 26.390573296152585}, {"from": "HSP90AB1", "to": "FGFR2", "weight": 25.649493574615366}, {"from": "HSP90AB1", "to": "GNAI2", "weight": 25.649493574615366}, {"from": "HSP90AB1", "to": "GNB2", "weight": 23.978952727983707}, {"from": "HSP90AB1", "to": "GNB3", "weight": 21.972245773362197}, {"from": "HSP90AB1", "to": "GRB10", "weight": 21.972245773362197}, {"from": "HSP90AB1", "to": "IGF2R", "weight": 23.02585092994046}, {"from": "HSP90AB1", "to": "IRAK1", "weight": 27.0805020110221}, {"from": "HSP90AB1", "to": "ITPR1", "weight": 21.972245773362197}, {"from": "HSP90AB1", "to": "JAK1", "weight": 23.978952727983707}, {"from": "HSP90AB1", "to": "KRAS", "weight": 24.849066497880003}, {"from": "HSP90AB1", "to": "KSR1", "weight": 21.972245773362197}, {"from": "HSP90AB1", "to": "LAMB2", "weight": 21.972245773362197}, {"from": "HSP90AB1", "to": "MAP2K7", "weight": 23.02585092994046}, {"from": "HSP90AB1", "to": "MAPK13", "weight": 23.978952727983707}, {"from": "HSP90AB1", "to": "MYO9B", "weight": 23.978952727983707}, {"from": "HSP90AB1", "to": "NCSTN", "weight": 24.849066497880003}, {"from": "HSP90AB1", "to": "PFKL", "weight": 23.978952727983707}, {"from": "HSP90AB1", "to": "PIM1", "weight": 21.972245773362197}, {"from": "HSP90AB1", "to": "PPARG", "weight": 24.849066497880003}, {"from": "HSP90AB1", "to": "PPP2R1A", "weight": 25.649493574615366}, {"from": "HSP90AB1", "to": "PTK2B", "weight": 23.02585092994046}, {"from": "HSP90AB1", "to": "RARA", "weight": 24.849066497880003}, {"from": "HSP90AB1", "to": "ROCK2", "weight": 24.849066497880003}, {"from": "HSP90AB1", "to": "SDC1", "weight": 23.02585092994046}, {"from": "HSP90AB1", "to": "STAT1", "weight": 25.649493574615366}, {"from": "HSP90AB1", "to": "TNK2", "weight": 23.02585092994046}, {"from": "HSP90AB1", "to": "TRIB3", "weight": 23.978952727983707}, {"from": "IRAK1", "to": "CHEK2", "weight": 27.725887222397812}, {"from": "IRAK1", "to": "DUSP10", "weight": 21.972245773362197}, {"from": "IRAK1", "to": "DUSP9", "weight": 21.972245773362197}, {"from": "IRAK1", "to": "ESR1", "weight": 26.390573296152585}, {"from": "IRAK1", "to": "FGFR2", "weight": 23.02585092994046}, {"from": "IRAK1", "to": "JAK1", "weight": 24.849066497880003}, {"from": "IRAK1", "to": "KALRN", "weight": 23.02585092994046}, {"from": "IRAK1", "to": "KRAS", "weight": 24.849066497880003}, {"from": "IRAK1", "to": "KSR1", "weight": 25.649493574615366}, {"from": "IRAK1", "to": "MAP2K7", "weight": 27.725887222397812}, {"from": "IRAK1", "to": "MAP4K2", "weight": 23.978952727983707}, {"from": "IRAK1", "to": "MAP4K3", "weight": 21.972245773362197}, {"from": "IRAK1", "to": "MAPK13", "weight": 28.903717578961647}, {"from": "IRAK1", "to": "MAST2", "weight": 21.972245773362197}, {"from": "IRAK1", "to": "PIM1", "weight": 26.390573296152585}, {"from": "IRAK1", "to": "PPARG", "weight": 25.649493574615366}, {"from": "IRAK1", "to": "PPP2R1A", "weight": 25.649493574615366}, {"from": "IRAK1", "to": "PTK2B", "weight": 25.649493574615366}, {"from": "IRAK1", "to": "PTPN23", "weight": 21.972245773362197}, {"from": "IRAK1", "to": "RARA", "weight": 23.978952727983707}, {"from": "IRAK1", "to": "ROCK2", "weight": 25.649493574615366}, {"from": "IRAK1", "to": "SNRK", "weight": 21.972245773362197}, {"from": "IRAK1", "to": "STAT1", "weight": 27.725887222397812}, {"from": "IRAK1", "to": "TLK2", "weight": 21.972245773362197}, {"from": "IRAK1", "to": "TNK2", "weight": 24.849066497880003}, {"from": "IRAK1", "to": "TRIB3", "weight": 25.649493574615366}, {"from": "ESR1", "to": "CHEK2", "weight": 25.649493574615366}, {"from": "ESR1", "to": "EGR2", "weight": 24.849066497880003}, {"from": "ESR1", "to": "FGFR2", "weight": 25.649493574615366}, {"from": "ESR1", "to": "GNAI2", "weight": 24.849066497880003}, {"from": "ESR1", "to": "GNB2", "weight": 21.972245773362197}, {"from": "ESR1", "to": "GRB10", "weight": 21.972245773362197}, {"from": "ESR1", "to": "KDM3A", "weight": 24.849066497880003}, {"from": "ESR1", "to": "KRAS", "weight": 25.649493574615366}, {"from": "ESR1", "to": "KSR1", "weight": 23.02585092994046}, {"from": "ESR1", "to": "MAP2K7", "weight": 21.972245773362197}, {"from": "ESR1", "to": "NCOA6", "weight": 21.972245773362197}, {"from": "ESR1", "to": "NCSTN", "weight": 23.02585092994046}, {"from": "ESR1", "to": "PEX14", "weight": 21.972245773362197}, {"from": "ESR1", "to": "PPARG", "weight": 32.58096538021482}, {"from": "ESR1", "to": "PPARGC1B", "weight": 21.972245773362197}, {"from": "ESR1", "to": "PPP2R1A", "weight": 26.390573296152585}, {"from": "ESR1", "to": "RARA", "weight": 31.780538303479457}, {"from": "ESR1", "to": "RNF14", "weight": 23.02585092994046}, {"from": "ESR1", "to": "ROCK2", "weight": 21.972245773362197}, {"from": "ESR1", "to": "RREB1", "weight": 23.02585092994046}, {"from": "ESR1", "to": "SMARCA2", "weight": 26.390573296152585}, {"from": "ESR1", "to": "SRF", "weight": 23.978952727983707}, {"from": "ESR1", "to": "STAT1", "weight": 32.18875824868201}, {"from": "ESR1", "to": "TNK2", "weight": 21.972245773362197}, {"from": "ESR1", "to": "TRIB3", "weight": 23.02585092994046}, {"from": "STAT1", "to": "CHEK2", "weight": 23.978952727983707}, {"from": "STAT1", "to": "EGR2", "weight": 23.978952727983707}, {"from": "STAT1", "to": "FGFR2", "weight": 23.02585092994046}, {"from": "STAT1", "to": "GNAI2", "weight": 23.02585092994046}, {"from": "STAT1", "to": "GNB2", "weight": 21.972245773362197}, {"from": "STAT1", "to": "GNB3", "weight": 23.02585092994046}, {"from": "STAT1", "to": "JAK1", "weight": 24.849066497880003}, {"from": "STAT1", "to": "KDM3A", "weight": 23.02585092994046}, {"from": "STAT1", "to": "KRAS", "weight": 27.725887222397812}, {"from": "STAT1", "to": "KSR1", "weight": 23.02585092994046}, {"from": "STAT1", "to": "MAP2K7", "weight": 25.649493574615366}, {"from": "STAT1", "to": "MAPK13", "weight": 21.972245773362197}, {"from": "STAT1", "to": "PIM1", "weight": 23.02585092994046}, {"from": "STAT1", "to": "PPARG", "weight": 30.91042453358316}, {"from": "STAT1", "to": "PPARGC1B", "weight": 21.972245773362197}, {"from": "STAT1", "to": "PPP2R1A", "weight": 28.33213344056216}, {"from": "STAT1", "to": "PTK2B", "weight": 27.0805020110221}, {"from": "STAT1", "to": "RARA", "weight": 30.91042453358316}, {"from": "STAT1", "to": "RNF14", "weight": 23.02585092994046}, {"from": "STAT1", "to": "RREB1", "weight": 21.972245773362197}, {"from": "STAT1", "to": "SMARCA2", "weight": 25.649493574615366}, {"from": "STAT1", "to": "SRF", "weight": 23.978952727983707}, {"from": "STAT1", "to": "TRIB3", "weight": 21.972245773362197}, {"from": "ROCK2", "to": "CAMK1G", "weight": 21.972245773362197}, {"from": "ROCK2", "to": "CDK18", "weight": 21.972245773362197}, {"from": "ROCK2", "to": "CDK19", "weight": 21.972245773362197}, {"from": "ROCK2", "to": "CHEK2", "weight": 24.849066497880003}, {"from": "ROCK2", "to": "GNAI2", "weight": 25.649493574615366}, {"from": "ROCK2", "to": "GNB2", "weight": 21.972245773362197}, {"from": "ROCK2", "to": "GNB3", "weight": 21.972245773362197}, {"from": "ROCK2", "to": "ITSN1", "weight": 23.02585092994046}, {"from": "ROCK2", "to": "JAK1", "weight": 23.02585092994046}, {"from": "ROCK2", "to": "KALRN", "weight": 25.649493574615366}, {"from": "ROCK2", "to": "KRAS", "weight": 23.02585092994046}, {"from": "ROCK2", "to": "KSR1", "weight": 23.978952727983707}, {"from": "ROCK2", "to": "MAP2K7", "weight": 21.972245773362197}, {"from": "ROCK2", "to": "MAPK13", "weight": 28.33213344056216}, {"from": "ROCK2", "to": "MAST2", "weight": 21.972245773362197}, {"from": "ROCK2", "to": "PIM1", "weight": 24.849066497880003}, {"from": "ROCK2", "to": "PPARG", "weight": 21.972245773362197}, {"from": "ROCK2", "to": "PTK2B", "weight": 23.978952727983707}, {"from": "ROCK2", "to": "SNRK", "weight": 21.972245773362197}, {"from": "ROCK2", "to": "TNK2", "weight": 23.978952727983707}, {"from": "ROCK2", "to": "TRIB3", "weight": 23.02585092994046}, {"from": "MAPK13", "to": "CDK18", "weight": 21.972245773362197}, {"from": "MAPK13", "to": "CDK19", "weight": 21.972245773362197}, {"from": "MAPK13", "to": "CHEK2", "weight": 23.02585092994046}, {"from": "MAPK13", "to": "GNAI2", "weight": 21.972245773362197}, {"from": "MAPK13", "to": "ITPR1", "weight": 23.02585092994046}, {"from": "MAPK13", "to": "JAK1", "weight": 24.849066497880003}, {"from": "MAPK13", "to": "KALRN", "weight": 23.978952727983707}, {"from": "MAPK13", "to": "KRAS", "weight": 23.978952727983707}, {"from": "MAPK13", "to": "KSR1", "weight": 25.649493574615366}, {"from": "MAPK13", "to": "MAP2K7", "weight": 27.725887222397812}, {"from": "MAPK13", "to": "MAP4K2", "weight": 23.978952727983707}, {"from": "MAPK13", "to": "MAP4K3", "weight": 23.02585092994046}, {"from": "MAPK13", "to": "MAST2", "weight": 21.972245773362197}, {"from": "MAPK13", "to": "MYO9B", "weight": 21.972245773362197}, {"from": "MAPK13", "to": "PIM1", "weight": 24.849066497880003}, {"from": "MAPK13", "to": "PTK2B", "weight": 26.390573296152585}, {"from": "MAPK13", "to": "SNRK", "weight": 23.02585092994046}, {"from": "MAPK13", "to": "TLK2", "weight": 23.02585092994046}, {"from": "MAPK13", "to": "TNK2", "weight": 21.972245773362197}, {"from": "MAPK13", "to": "TRIB3", "weight": 23.02585092994046}, {"from": "PPARG", "to": "CHEK2", "weight": 23.978952727983707}, {"from": "PPARG", "to": "EGR2", "weight": 26.390573296152585}, {"from": "PPARG", "to": "GNAI2", "weight": 23.978952727983707}, {"from": "PPARG", "to": "GNB2", "weight": 23.02585092994046}, {"from": "PPARG", "to": "GNB3", "weight": 21.972245773362197}, {"from": "PPARG", "to": "KDM3A", "weight": 23.02585092994046}, {"from": "PPARG", "to": "KRAS", "weight": 23.02585092994046}, {"from": "PPARG", "to": "MAP2K7", "weight": 23.02585092994046}, {"from": "PPARG", "to": "NCOA6", "weight": 25.649493574615366}, {"from": "PPARG", "to": "PPARGC1B", "weight": 24.849066497880003}, {"from": "PPARG", "to": "PPP2R1A", "weight": 24.849066497880003}, {"from": "PPARG", "to": "PTK2B", "weight": 21.972245773362197}, {"from": "PPARG", "to": "RARA", "weight": 32.18875824868201}, {"from": "PPARG", "to": "RNF14", "weight": 23.02585092994046}, {"from": "PPARG", "to": "RREB1", "weight": 23.02585092994046}, {"from": "PPARG", "to": "SMARCA2", "weight": 27.725887222397812}, {"from": "PPARG", "to": "SRF", "weight": 23.978952727983707}, {"from": "PPARG", "to": "TNK2", "weight": 21.972245773362197}, {"from": "PPARG", "to": "TRIB3", "weight": 23.02585092994046}, {"from": "PTK2B", "to": "FGFR2", "weight": 23.02585092994046}, {"from": "PTK2B", "to": "GNAI2", "weight": 26.390573296152585}, {"from": "PTK2B", "to": "GNB2", "weight": 24.849066497880003}, {"from": "PTK2B", "to": "GNB3", "weight": 23.978952727983707}, {"from": "PTK2B", "to": "ITPR1", "weight": 23.978952727983707}, {"from": "PTK2B", "to": "JAK1", "weight": 27.725887222397812}, {"from": "PTK2B", "to": "KRAS", "weight": 24.849066497880003}, {"from": "PTK2B", "to": "KSR1", "weight": 23.02585092994046}, {"from": "PTK2B", "to": "MAP2K7", "weight": 27.725887222397812}, {"from": "PTK2B", "to": "PIM1", "weight": 23.02585092994046}, {"from": "PTK2B", "to": "PPP2R1A", "weight": 23.978952727983707}, {"from": "PTK2B", "to": "RARA", "weight": 23.978952727983707}, {"from": "PTK2B", "to": "TNK2", "weight": 25.649493574615366}, {"from": "PTK2B", "to": "TRIB3", "weight": 23.02585092994046}, {"from": "PPP2R1A", "to": "GNAI2", "weight": 23.978952727983707}, {"from": "PPP2R1A", "to": "GNB2", "weight": 23.978952727983707}, {"from": "PPP2R1A", "to": "GNB3", "weight": 23.978952727983707}, {"from": "PPP2R1A", "to": "ITPR1", "weight": 23.978952727983707}, {"from": "PPP2R1A", "to": "JAK1", "weight": 25.649493574615366}, {"from": "PPP2R1A", "to": "KRAS", "weight": 23.978952727983707}, {"from": "PPP2R1A", "to": "KSR1", "weight": 23.02585092994046}, {"from": "PPP2R1A", "to": "LAMA5", "weight": 21.972245773362197}, {"from": "PPP2R1A", "to": "MAP2K7", "weight": 23.978952727983707}, {"from": "PPP2R1A", "to": "NCSTN", "weight": 23.02585092994046}, {"from": "PPP2R1A", "to": "RARA", "weight": 23.02585092994046}, {"from": "PPP2R1A", "to": "SDC1", "weight": 24.849066497880003}, {"from": "PPP2R1A", "to": "TRIB3", "weight": 21.972245773362197}, {"from": "CHEK2", "to": "CDK18", "weight": 23.02585092994046}, {"from": "CHEK2", "to": "FGFR2", "weight": 21.972245773362197}, {"from": "CHEK2", "to": "KSR1", "weight": 21.972245773362197}, {"from": "CHEK2", "to": "MAP2K7", "weight": 21.972245773362197}, {"from": "CHEK2", "to": "PIM1", "weight": 23.978952727983707}, {"from": "CHEK2", "to": "RARA", "weight": 23.02585092994046}, {"from": "CHEK2", "to": "SMARCA2", "weight": 21.972245773362197}, {"from": "CHEK2", "to": "SNRK", "weight": 21.972245773362197}, {"from": "CHEK2", "to": "TLK2", "weight": 21.972245773362197}, {"from": "CHEK2", "to": "TNK2", "weight": 23.978952727983707}, {"from": "CHEK2", "to": "TRIB3", "weight": 21.972245773362197}, {"from": "RARA", "to": "EGR2", "weight": 23.978952727983707}, {"from": "RARA", "to": "GNAI2", "weight": 23.02585092994046}, {"from": "RARA", "to": "GNB2", "weight": 21.972245773362197}, {"from": "RARA", "to": "GNB3", "weight": 21.972245773362197}, {"from": "RARA", "to": "KDM3A", "weight": 23.02585092994046}, {"from": "RARA", "to": "NCOA6", "weight": 23.02585092994046}, {"from": "RARA", "to": "PPARGC1B", "weight": 23.02585092994046}, {"from": "RARA", "to": "SMARCA2", "weight": 27.0805020110221}, {"from": "RARA", "to": "SRF", "weight": 23.978952727983707}, {"from": "KSR1", "to": "JAK1", "weight": 24.849066497880003}, {"from": "KSR1", "to": "KALRN", "weight": 23.978952727983707}, {"from": "KSR1", "to": "MAP2K7", "weight": 24.849066497880003}, {"from": "KSR1", "to": "MAP4K2", "weight": 21.972245773362197}, {"from": "KSR1", "to": "PIM1", "weight": 23.02585092994046}, {"from": "KSR1", "to": "TNK2", "weight": 23.978952727983707}, {"from": "KSR1", "to": "TRIB3", "weight": 21.972245773362197}, {"from": "MAP2K7", "to": "FGFR2", "weight": 23.02585092994046}, {"from": "MAP2K7", "to": "JAK1", "weight": 25.649493574615366}, {"from": "MAP2K7", "to": "KRAS", "weight": 23.02585092994046}, {"from": "MAP2K7", "to": "PIM1", "weight": 25.649493574615366}, {"from": "MAP2K7", "to": "TNK2", "weight": 23.02585092994046}, {"from": "GNAI2", "to": "GNB2", "weight": 30.91042453358316}, {"from": "GNAI2", "to": "GNB3", "weight": 30.91042453358316}, {"from": "GNAI2", "to": "ITPR1", "weight": 26.390573296152585}, {"from": "GNAI2", "to": "ITSN1", "weight": 21.972245773362197}, {"from": "GNAI2", "to": "KALRN", "weight": 23.02585092994046}, {"from": "GNAI2", "to": "KRAS", "weight": 26.390573296152585}, {"from": "GNAI2", "to": "SPR", "weight": 21.972245773362197}, {"from": "KRAS", "to": "FGFR2", "weight": 23.02585092994046}, {"from": "KRAS", "to": "GNB2", "weight": 25.649493574615366}, {"from": "KRAS", "to": "GNB3", "weight": 23.978952727983707}, {"from": "KRAS", "to": "ITPR1", "weight": 23.02585092994046}, {"from": "KRAS", "to": "JAK1", "weight": 23.978952727983707}, {"from": "JAK1", "to": "FGFR2", "weight": 23.02585092994046}, {"from": "JAK1", "to": "GNB2", "weight": 23.02585092994046}, {"from": "JAK1", "to": "PIM1", "weight": 23.02585092994046}, {"from": "JAK1", "to": "TNK2", "weight": 24.849066497880003}, {"from": "JAK1", "to": "TRIB3", "weight": 21.972245773362197}, {"from": "PIM1", "to": "MAST2", "weight": 21.972245773362197}, {"from": "PIM1", "to": "TLK2", "weight": 21.972245773362197}, {"from": "PIM1", "to": "TNK2", "weight": 21.972245773362197}, {"from": "PIM1", "to": "TRIB3", "weight": 23.02585092994046}, {"from": "TNK2", "to": "FGFR2", "weight": 23.978952727983707}, {"from": "GNB2", "to": "GNB3", "weight": 31.354942159291497}, {"from": "GNB2", "to": "ITPR1", "weight": 26.390573296152585}, {"from": "GNB3", "to": "ITPR1", "weight": 25.649493574615366}, {"from": "ITPR1", "to": "IGF2R", "weight": 21.972245773362197}, {"from": "ITPR1", "to": "NCSTN", "weight": 21.972245773362197}, {"from": "ITPR1", "to": "SDC1", "weight": 21.972245773362197}, {"from": "NCSTN", "to": "IGF2R", "weight": 24.849066497880003}, {"from": "NCSTN", "to": "ITGA1", "weight": 21.972245773362197}, {"from": "NCSTN", "to": "SDC1", "weight": 21.972245773362197}, {"from": "SMARCA2", "to": "EGR2", "weight": 23.02585092994046}, {"from": "SMARCA2", "to": "PPARGC1B", "weight": 21.972245773362197}, {"from": "EGR2", "to": "KDM3A", "weight": 21.972245773362197}, {"from": "EGR2", "to": "NCOA6", "weight": 21.972245773362197}, {"from": "EGR2", "to": "SRF", "weight": 21.972245773362197}, {"from": "KALRN", "to": "ITSN1", "weight": 23.02585092994046}, {"from": "KDM3A", "to": "SRF", "weight": 23.02585092994046}, {"from": "IGF2R", "to": "CD46", "weight": 21.972245773362197}, {"from": "IGF2R", "to": "SERPINA1", "weight": 21.972245773362197}, {"from": "LAMB2", "to": "COL6A3", "weight": 26.390573296152585}, {"from": "LAMB2", "to": "ITGA1", "weight": 21.972245773362197}, {"from": "LAMB2", "to": "LAMA5", "weight": 23.978952727983707}, {"from": "MAP4K2", "to": "MAP4K3", "weight": 23.02585092994046}, {"from": "ITGA1", "to": "ITGA7", "weight": 23.978952727983707}, {"from": "LAMA5", "to": "COL6A3", "weight": 23.02585092994046}, {"from": "DUSP9", "to": "DUSP10", "weight": 23.978952727983707}, {"from": "SLCO1A2", "to": "SLCO1B3", "weight": 21.972245773362197}]);

        // adding nodes and edges to the graph
        data = {nodes: nodes, edges: edges};

        var options = {
    "configure": {
        "enabled": true,
        "filter": [
            "nodes",
            "physics"
        ]
    },
    "edges": {
        "color": {
            "inherit": true
        },
        "smooth": {
            "enabled": false,
            "type": "continuous"
        }
    },
    "interaction": {
        "dragNodes": true,
        "hideEdgesOnDrag": false,
        "hideNodesOnDrag": false
    },
    "physics": {
        "enabled": true,
        "repulsion": {
            "centralGravity": 0.2,
            "damping": 0.09,
            "nodeDistance": 120,
            "springConstant": 0.05,
            "springLength": 220
        },
        "solver": "repulsion",
        "stabilization": {
            "enabled": true,
            "fit": true,
            "iterations": 1000,
            "onlyDynamicEdges": false,
            "updateInterval": 50
        }
    }
};
        
        

        
        // if this network requires displaying the configure window,
        // put it in its div
        options.configure["container"] = document.getElementById("config");
        

        network = new vis.Network(container, data, options);
	 
        


        

        return network;

    }

    drawGraph();

</script>
</body>
</html>

```