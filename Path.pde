/**
 * Basic Path class. A path contains a collection of Interaction elements.
 * Path class is used for querying shortest path between any two user
 * provided proteins.
 *
 * A simple (perhaps inefficient) algorithm is used for shortest path.
 *
 * @author Pranav Srinivas
 * @version 23 May 2015
 */

class Path
{
  ArrayList<Interaction> elements; // Interaction elements of the path
  
  Path()
  {
    elements = new ArrayList<Interaction>();
  }
  
  // Construct from other path
  Path(Path other)
  {
    elements = new ArrayList<Interaction>();
    for (Interaction link : other.elements)
    {
      elements.add(link);
    }
  }
  
  // Get Path length
  int getLength()
  {
    return elements.size();
  }
  
  void clear()
  {
    elements.clear();
  }
  
  void setHighlight(boolean value)
  {
    for (Interaction link : elements)
      link.setHighlight(value);
  }
  
  void add(Interaction link)
  {
    elements.add(link);
  }
  
  /** 
   * Remove last interaction on the path
   */
  void removeLast()
  {
    if (elements.size() > 0) 
    {
      elements.remove(elements.size() - 1);
    }
  }
}

/**
 * Determine shortest path between two proteins (user provided)
 */
void determineShortestPath()
{
  if (menuID != 8) return;
  
  // Nothing to do if user didn't type anything
  if (savedUserText.equals("")) return;
  
  String [] prots = savedUserText.split(" +");  // Split by one or more spaces
  
  // Shortest path can be queried only between two proteins
  if (prots.length != 2) return;
  
  Protein startProtein  = nameToProteinMap.get( prots[0] );  // Lookup by name
  Protein targetProtein = nameToProteinMap.get( prots[1] );  // Lookup by name
  
  // Return if start or end protein can't be found
  if (startProtein == null || targetProtein == null) return;
  
  // Clear all existing highlights  
  for (Interaction link : interactions)
  {
    link.setHighlight(false);
  }
  
  Path currentPath = new Path();
  ArrayList<Path> allPaths = new ArrayList<Path>();
  
  // Find all paths between startProtein and targetProtein
  // and look for the shortest one. There are better algorithms
  // for shortest path
  queryPaths(startProtein, targetProtein, currentPath, allPaths);

  // Pick the shortest from all paths (naive approach)
  Path shortestPath = null;
  if (allPaths.size() > 0) 
  {
    shortestPath = allPaths.get(0);  
    for (Path path : allPaths)
    {
      if (path.getLength() < shortestPath.getLength())
        shortestPath = path;
    }
  }
  
  // Highlight shortest path
  if (shortestPath != null) shortestPath.setHighlight(true);
  
  allPaths.clear();
  
  // Clear all isVisited flag
  for (Protein p : proteins)
  {
    p.setIsVisited(false);
  }
  for (Interaction link : interactions)
  {
    link.setIsVisited(false);
  }
}

/**
 * Recursively find all paths between 'from' and 'target'.
 * Use 'isVisited' flag  to not go over same link more than once
 */
void queryPaths(Protein from, Protein target, Path currentPath, ArrayList<Path> allPaths)
{
  // Reached target! Add currentPath to path list and return
  if (from == target) 
  {
    Path pth = new Path(currentPath);
    allPaths.add(pth);
    return;
  }           
   
  // Go through each interaction of 'from' protein  
  for (Interaction link : from.getLinks())
  {
    Protein other = (from == link.getSourceProtein()) ? link.getSinkProtein() : link.getSourceProtein();
    if (link.getIsVisited() == true) continue; // Already passed through 
    link.setIsVisited(true);
    currentPath.add(link);
    queryPaths(other, target, currentPath, allPaths);  // recursive call
    currentPath.removeLast();
    link.setIsVisited(false);
  }
}
