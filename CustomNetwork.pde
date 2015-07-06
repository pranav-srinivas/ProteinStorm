/**
 * This file contains code for visualizing custom biological networks.
 * The networks are in ASCII format (two proteins separated by tab):
 * For Example:
 *                ProteinA  ProteinB 
 *                ProteinC  ProteinD
 * Also contains a method for creating a random network of 100 proteins.
 *
 * @author Pranav Srinivas
 * @version 18 May 2015
 */


String [] customFiles = {"CN1.txt", "CN2.txt", "CN3.txt"}; // custom network file names
int cIndex = 0;

/** 
 * Visualize custom biological network (tab separated protein names)
 */
void loadCustomBiologicalNetwork()
{
  clearNetwork(); // Clear any existing network;
  
  // Load cutsom network file as array of strings
  String [] lines = loadStrings(customFiles[cIndex]);

  float radius = 9;
  
  // Process each line creating proteins and interactions
  for (int i = 0; i < lines.length; i++) 
  {
    String [] prots = split( lines[i], '\t');  // Split by tab
    if (prots.length == 2)
    {
      String name1 = prots[0];  // first protein
      String name2 = prots[1];  // second protein
      
      Protein protein1 = nameToProteinMap.get(name1); // lookup
      if (protein1 == null)
      {
        protein1 = new Protein(name1);
        nameToProteinMap.put(name1, protein1);
        protein1.setRadius(100);
        protein1.setRepulsionStrength(-5);
        protein1.setPositionLimit(radius, radius, width-radius, height-radius);
        proteins.add(protein1);
      }
      
      Protein protein2 = nameToProteinMap.get(name2);  // lookup
      if (protein2 == null)
      {
        protein2 = new Protein(name2);
        nameToProteinMap.put(name2, protein2);
        protein2.setRadius(100);
        protein2.setRepulsionStrength(-5);
        protein2.setPositionLimit(radius, radius, width-radius, height-radius);
        proteins.add(protein2);      
      }

      Interaction link = new Interaction(protein1, protein2);
      link.setLength(20);
      link.setRigidity(1);
      interactions.add(link);  
      protein1.addLink(link);
      protein2.addLink(link);
    }
  }
  
  cIndex++;             // next custom network file for future selection 
  cIndex = cIndex % 3;  // Recycle!, only 3 available for now
}


/** 
 * Visualize Random biological network
 */
void loadRandomBiologicalNetwork()
{
  clearNetwork(); // Clear any existing network;
  
  // Randomly create proteins
  float radius = proteinDiameter/2;
  for (int i = 0; i < 100; i++) 
  {
    Protein protein = new Protein(width/2+random(-200, 200), height/2+random(-200, 200));
    protein.setPositionLimit(radius, radius, width-radius, height-radius);
    protein.setRadius(100);
    protein.setRepulsionStrength(-12);
    proteins.add(protein);
  }

  // Randomly create interactions
  for (int j = 0 ; j < proteins.size()-1; j++) 
  {
    int count = floor(random(1, 3));
    for (int i = 0; i < count; i++) 
    {
      int r = floor(random(j+1, proteins.size()));
      Interaction link = new Interaction(proteins.get(j), proteins.get(r));
      link.setLength(20);
      link.setRigidity(1);
      interactions.add(link);
    }
  }
}
