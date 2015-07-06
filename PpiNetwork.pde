/**
 * This file contains Protein Protein Interaction Network class.
 * One background Protein-Protein Interaction (PPI) Network is
 * created from downloaded data in JSON format (source: HPRD).
 * The background PPI network is huge (19653 proteins, 39240 interactions)
 * and is never displayed.
 *
 * The background PPI network is used for enriching a user provided gene list.
 * The enriched network (subnetwork of background PPI network) is much smaller
 * and is used for visualization.
 *
 * @author Pranav Srinivas
 * @version 17 May 2015
 */

class PpiNetwork
{
  private ArrayList<Protein>       PpiProteins;     // Proteins of the PPI network
  private ArrayList<Interaction>   PpiInteractions; // Interactions of the PPI network
  private HashMap<String, Protein> idToProtein;     // Dictionary map from id to protein
  private HashMap<String, Protein> nameToProtein;   // Dictionary map from name to protein

  // Constructor
  PpiNetwork()
  {
    PpiProteins     = new ArrayList<Protein>();
    PpiInteractions = new ArrayList<Interaction>();
    idToProtein     = new HashMap<String, Protein>();    
    nameToProtein   = new HashMap<String, Protein>();
  }
  
  // Constructor
  PpiNetwork(JSONArray Jnodes, JSONArray Jlinks)
  {
    PpiProteins     = new ArrayList<Protein>();
    PpiInteractions = new ArrayList<Interaction>();
    idToProtein     = new HashMap<String, Protein>();
    nameToProtein   = new HashMap<String, Protein>();
    
    createProteins(Jnodes);      // Create all protein nodes
    createInteractions(Jlinks);  // Create all interactions
  }
  
   
  //-------------------------------------------------------------
  // Accessor Methods (Getters)
  //-------------------------------------------------------------
  
  // Get total protein count
  int getProteinCount()
  {
    return PpiProteins.size();
  }
  
  // Get total interaction count
  int getInteractionCount()
  {
    return PpiInteractions.size();
  }  
  
  /**
    * Get protein with a given Id
    * @param id is the string id of the protein
    * @return Protein object corresponding to the id
    */
  Protein getProtein(String id)
  {
    return idToProtein.get(id);
  }
  
  /**
    * Get protein with a given name
    * @param pname is the name of the protein
    * @return Protein object corresponding to the name 'pname'
    */
  Protein getProteinFromName(String pname)
  {
    return nameToProtein.get(pname);
  }
  
  // Get protein list
  ArrayList<Protein> getProteinList()
  {
    return PpiProteins;
  }
  
  // Get interaction list
  ArrayList<Interaction> getInteractionList()
  {
    return PpiInteractions;
  }
  
  /**
   * Create all PPI proteins from the passed JSON Array
   */
  void createProteins(JSONArray Jnodes)
  {
    float radius = proteinDiameter/2;
    for (int i = 0; i < Jnodes.size(); i++)
    {
      JSONObject Jnode = Jnodes.getJSONObject(i);
      String pname     = Jnode.getString("geneSymbol");
      String pid       = Jnode.getString("hprd_id");
    
      // Create unique protein object and set all attributes
      Protein protein = new Protein(width/2+random(-200, 200), height/2+random(-200, 200));
      protein.setPositionLimit(radius, radius, width-radius, height-radius);
      protein.setRadius(100);
      protein.setRepulsionStrength(-10);
      protein.setName(pname);
      protein.setId(pid);
     
      idToProtein.put(pid, protein);      // Add protein to dictionary for easy retrieval
      nameToProtein.put(pname, protein);  // Add protein to dictionary for easy retrieval
      PpiProteins.add(protein);           // Add to ArrayList of proteins 
    }
  }
  
  
  /**
   * Create all PPI interactions from the passed JSON Array
   */
  void createInteractions(JSONArray Jlinks)
  {
    for (int i = 0; i < Jlinks.size(); i++)
    {
      JSONObject Jlink = Jlinks.getJSONObject(i);
      String id1       = Jlink.getString("interactor_1_hprd_id");
      String id2       = Jlink.getString("interactor_2_hprd_id");
      
      Protein proteinOne = getProtein( id1 );
      Protein proteinTwo = getProtein( id2 );
      Interaction link = new Interaction(proteinOne, proteinTwo);
      link.setLength(25);
      link.setRigidity(0.4);
      
      PpiInteractions.add(link);  // Add to ArrayList of interactions   
    }
  }
  
  /**
   * Clear the PPI network
   */
  void clear()
  {
    PpiProteins.clear();
    PpiInteractions.clear();
    idToProtein.clear();
    nameToProtein.clear();
  }
}


PpiNetwork backgroundPpiNetwork;  // global Background HPRD PPI Network

/**
 * Load full Background HPRD Protein-Protein Interaction (PPI) Network
 * The PPI network is represented as undirected graph.
 * The PPI network is obtained from HPRD website in XML format and later
 * converted into JSON format using xml2json utility. The full PPI network
 * has 19653 proteins (nodes) and 39240 interactions (edges)
 * Networks are enriched from the Background PPI network based on a gene set
 */
void loadBackgroundPPINetwork()
{
  if (backgroundPpiNetwork != null) return;  // Already loaded
  
  // Load Background PPI JSON file
  JSONObject ppis = loadJSONObject("BackgroundPPI.json");
  
  // Get all protein nodes from Background PPI JSON
  JSONArray Jnodes = ppis.getJSONArray("Nodes");
  
  // Get all interaction links from Background PPI JSON
  JSONArray Jlinks = ppis.getJSONArray("Links");
  
  // Create PPI network
  backgroundPpiNetwork = new PpiNetwork(Jnodes, Jlinks);
}

String [] geneSetFiles = {"GS1.txt", "GS2.txt", "GS3.txt"};
int gsIndex = 0;

/**
 * Derive sub-network (enrich) from the background PPI network based 
 * on user provided gene set.
 */
void enrichPPIFromGeneSet()
{
  clearNetwork(); // Clear any existing displayed network;
  
  // Load Background PPI network if not done yet
  if (backgroundPpiNetwork == null)
  {
    loadBackgroundPPINetwork();
  }
  
  // Load Gene Set file as array of strings
  String [] lines = loadStrings(geneSetFiles[gsIndex]);  
  
  // Create and populate geneList from the file
  ArrayList<String> geneList = new ArrayList<String>();
  for (int i = 0; i < lines.length; i++) 
  {
    String [] genes = split( lines[i], '\t');  // Split by tab
    for (int j = 0; j < genes.length; j++)
    {
      geneList.add( genes[j] );
    }
  }
  
  enrichPPIFromGeneList(geneList);
  
  gsIndex++;              // next gene set file for future selection 
  gsIndex = gsIndex % 3;  // Recycle!, only 3 available for now
}
   
   
/**
 * Helper method for enrichPPIFromGeneSet() and enrichPPIFromGeneExpression().
 * Derives sub-network from the background PPI network based on array list
 * of gene names.
 */   
void enrichPPIFromGeneList(ArrayList<String> geneList)
{ 
  // Array list for enriched proteins and interactions
  ArrayList<Protein>     enrichedProteins     = new ArrayList<Protein>();
  ArrayList<Interaction> enrichedInteractions = new ArrayList<Interaction>();
 
  // Query relevant proteins from the background PPI network and add to
  // enriched proteins collection. Also mark protein as enriched.
  for (String geneName : geneList)
  {
    Protein p = backgroundPpiNetwork.getProteinFromName(geneName);
    if (p != null) 
    {
      p.setIsEnriched(true);
      enrichedProteins.add(p);
    }
  }
  
  // Add all interactions from the background PPI network whose source or sink 
  // proteins are marked as enriched, to enriched interactions collection.
  for (Interaction link : backgroundPpiNetwork.getInteractionList())
  {
    Protein source = link.getSourceProtein();
    Protein sink   = link.getSinkProtein();
    if (source.getIsEnriched() == true || sink.getIsEnriched() == true)
    {
      enrichedInteractions.add(link);     
      source.incrementLinksCount();
      sink.incrementLinksCount();
    }
  }
  
  // Further refine enriched interactions based on link count
  for (int i = 0; i < enrichedInteractions.size(); i++)
  {
    Interaction link = enrichedInteractions.get(i);
    Protein source   = link.getSourceProtein();
    Protein sink     = link.getSinkProtein();  
    
    if (source.getLinksCount() < 2 || sink.getLinksCount() < 2) 
    {
      if (!source.getIsEnriched()) source.clearLinksCount();
      if (!sink.getIsEnriched())   sink.clearLinksCount();
      enrichedInteractions.remove(i);
      i--;
      continue;
    }
    
    if (!source.getIsEnriched())
    {
      source.setIsEnriched(true);
      enrichedProteins.add(source);
    }
    if (!sink.getIsEnriched())
    {
      sink.setIsEnriched(true);
      enrichedProteins.add(sink);
    }
  }

  // Add enriched proteins to map for lookup by name  
  for (Protein p : enrichedProteins)
  {
    nameToProteinMap.put(p.getName(), p);
    p.reset();
  }
  
  // Populate links for proteins
  for (Interaction link : enrichedInteractions)
  {
    link.getSourceProtein().addLink(link);
    link.getSinkProtein().addLink(link);
  }
  
  // Assign to drawn proteins and interactions list
  proteins     = enrichedProteins;
  interactions = enrichedInteractions;
}

String [] geneExpressionFiles = {"GE1.txt", "GE2.txt"};
int geIndex = 0;
      
/**
 * Derive sub-network from the background PPI network based on genes in the
 * gene expression data. Also perform animation
 */      
void enrichPPIFromGeneExpression()
{
  clearNetwork(); // Clear any existing displayed network;
  
  // Load Background PPI network if not done yet
  if (backgroundPpiNetwork == null)
  {
    loadBackgroundPPINetwork();
  }
  
  // Load Gene Expression file as array of strings
  String [] lines = loadStrings( geneExpressionFiles[geIndex] );
  
  // Create and populate geneList from the file
  ArrayList<String> geneList = new ArrayList<String>();
  for (int i = 0; i < lines.length; i++) 
  {
    String [] genes = split( lines[i], '\t');
    for (int j = 0; j < genes.length; j++)
    {
      geneList.add( genes[j] );
    }
    break;
  }
  
  enrichPPIFromGeneList(geneList);
  
  geIndex++;              // next gene expression file for future selection 
  geIndex = geIndex % 2;  // Recycle!, only 2 available for now
}

