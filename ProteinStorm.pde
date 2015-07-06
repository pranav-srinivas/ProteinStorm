/**
 * This is the top level main program file for ProteinStorm.
 * ProteinStorm uses force layout for visualization of
 * Biological networks.
 *
 * Protein Storm features include
 * (a) Visualization of user created custom networks
 * (b) Visulization of enriched Protein-Protein Interaction (PPI) networks
 * (c) Visualization of Gene Expression data
 * (d) Identification and Visualization of Hub proteins
 * (e) Shortest Path identification and visualization
 * 
 * The program primarily uses two classes: 'Protein' and 'Interaction'
 * A network file is simply a series of tab separated protein names in ASCII format
 * 
 * Background PPI network was obtained from Human Protein Reference Database (HPRD) (http://hprd.org)
 *
 * @author Pranav Srinivas
 * @version 14 May 2015
 */

import java.util.ArrayList;

ArrayList<Protein>         proteins = new ArrayList<Protein>();             // Visualized proteins
ArrayList<Interaction> interactions = new ArrayList<Interaction>();         // Visualized interactions
HashMap<String, Protein> nameToProteinMap = new HashMap<String, Protein>(); // For efficient lookup

float proteinDiameter = 18;              // Diameter in pixel of visualized protein
float zoomLevel = 1.0;                   // Zoom level control
int buttonWidth  = 100;                  // Width  of menu buttons (in pixels)
int buttonHeight = 40;                   // Height of menu buttons (in pixels)
int menuID       = 1;                    // Menu ID number (1,2....8)
float  prevMouseX = 0,  prevMouseY = 0;  // Previous mouse location for dragging
float deltaMouseX = 0, deltaMouseY = 0;  // Delta mouse dragging
float panX = 0, panY = 0;                // Pan Control
static final int HUB_THRESHOLD = 8;      // Threshold for hub proteins
PFont tFont;                             // Font for drawing texts
String typing        = "";               // text currently being typed
String savedUserText = "";               // saved text when return is hit

/**
 * Processing display setup
 */

void setup()
{
  size(1400, 800);
  background(255);
  smooth();
  noStroke();
  
  colorMode(RGB, 255, 255, 255);
  tFont = createFont("Arial",12,true);
  textFont(tFont);
  
  loadRandomBiologicalNetwork();
}

/**
 * Processing display draw
 */

void draw()
{
  background(0);
  displayMenuButtons();  // Menu byttons at the top left
  
  if (menuID == 8) displayUserTextInputArea();  // Text for Shortest Path

  // Zoom to desired level
  translate(width/2 + panX, height/2 + panY);
  scale(zoomLevel);
  translate(-width/2 + panX, -height/2 + panY);
  
  // Visualize network
  viewNetwork();
}

/**
 * Display user typed texts
 */
void displayUserTextInputArea()
{  
  String s = savedUserText.equals("") ? typing : savedUserText;
  if (s.equals("")) return;
  
  fill(0, 0, 0);
  textFont(tFont, 14);
  textAlign(LEFT);
  text(s, 8*buttonWidth + 20, 20);
}

/**
 * Draw rectangles and texts for menu buttons
 */

void displayMenuButtons()
{
  stroke(0);
  fill (255, 215, 0);
  for (int i = 0; i < 8; i++) 
  {
    if (i == menuID - 1) fill(174, 221, 60);
    else                 fill(148, 73, 232);
    rect(i*buttonWidth, 0, buttonWidth, buttonHeight);
  }
  
  fill(132, 112, 255);
  rect(8*buttonWidth, 0, 2*buttonWidth, buttonHeight);
  
  fill(0, 0, 0);
  textFont(tFont, 14);
  textAlign(CENTER);

  String[] menuItems = { "Random",       "Network", 
                         "Custom",       "Network",
                         "PPI Enriched", "Network",
                         "Gene",         "Expression",
                         "Hub",          "Proteins",
                         "Fit",          "View", 
                         "Protein",      "Structure",
                         "Shortest",     "Path"
                       };
  
  for (int i = 0; i < 16; i++) 
  {
    int xPos = (i/2)*buttonWidth + buttonWidth/2;
    int yPos = (i % 2 == 0) ? 15 : 30;
    text(menuItems[i], xPos, yPos);
  }
}


/**
 * Clear all proteins and interactions
 */

void clearNetwork()
{
  // Clear various flags on proteins
  for (Protein p : proteins)
  {
    p.setIsEnriched(false);
    p.clearLinksCount();
    p.setIsHub(false);
    p.setIsVisited(false);
    p.clearLinks();
  }
  
  // Clear various flags on interactions
  for (Interaction link : interactions)
  {
    link.setHighlight(false);
  }
  
  // Clear existing proteins, interactions and the map
  proteins.clear();
  interactions.clear();
  nameToProteinMap.clear();
  
  savedUserText = "";
  typing        = ""; 
}

/**
 * Visualize network based on user menu selection
 */

void viewNetwork()
{ 
  // Proteins repulse each other!
  for (Protein protein : proteins)
  {
    protein.repulse(proteins);
  }

  // All protein interactions apply their forces
  for (Interaction interaction : interactions)
  {
    interaction.applyForces();
  }

  // Update velocity and position vectors for proteins
  for (Protein protein : proteins)
  {
    protein.updateVelocityAndPosition();
  }

  // Draw Interactions
  for (Interaction interaction : interactions)
  {
    if (interaction.getHighlight())
    {      
      stroke(255, 165, 0);
      strokeWeight(2);
    }
    else
    {
      stroke(0, 130, 164);
      strokeWeight(1);
    }
    line(interaction.getSourceProtein().getX(), interaction.getSourceProtein().getY(),
         interaction.getSinkProtein().getX(),   interaction.getSinkProtein().getY());
  }

  // Draw Proteins
  noStroke();
  textFont(tFont, 4);
  textAlign(CENTER);
 
  for (int i = 0; i < proteins.size(); i++)
  {
    Protein protein = proteins.get(i);

    int clr = i % 10;
    switch( clr ) {
      case 0: fill(255, 69, 0); break;
      case 1: fill(127, 255, 0); break;
      case 2: fill(132, 112, 255); break;
      case 3: fill(246, 255, 0); break;
      case 4: fill(255, 215, 0); break;
      case 5: fill(205, 92, 92); break;
      case 6: fill(255, 127, 80); break;
      case 7: fill(218, 112, 214); break;
      case 8: fill(127, 255, 212); break;
      case 9: fill(0, 206, 209); break;
      default: fill(100, 100, 100); break;
    }
    
    if (protein.getIsHub())
    {
      ellipse(protein.getX(), protein.getY(), 1.4*proteinDiameter, 1.4*proteinDiameter);
    }
    else
    {
      ellipse(protein.getX(), protein.getY(), proteinDiameter-4, proteinDiameter-4);
    }

    if (!protein.getName().equals(""))
    {
      fill(0);
      text(protein.getName(), protein.getX(), protein.getY());
    }
  }  
}


/**
 * Mark Hub proteins based on interaction count in the enriched network
 */

void determineHubProteins()
{ 
  for (Protein p : proteins)
  {
    if (p.getLinksCount() > HUB_THRESHOLD)
    {
      p.setIsHub(true);
    }
  } 
}


/**
 * Mouse events (mouseWheel for zoom in/out)
 * Control zoom level for Zoomin and ZoomOut
 */ 
void mouseWheel(MouseEvent evt) 
{
  float eventCount = evt.getCount();
  if (eventCount > 0) zoomLevel += 0.05;
  else                zoomLevel -= 0.05;  
  zoomLevel = constrain(zoomLevel, 0.2, 12.0);
}


/** 
 * Track drag direction
 */
void mouseDragged() 
{
  deltaMouseX = mouseX - prevMouseX;
  deltaMouseY = mouseY - prevMouseY;
}

/** 
 * Pan based on drag direction
 */
void mouseReleased() 
{
  panX += deltaMouseX;
  panY += deltaMouseY;  
  deltaMouseX = 0;
  deltaMouseY = 0;
}

/** 
 * Menu button selection and handling
 */
void mousePressed()
{
  int menuIndex = mouseX / buttonWidth;

  // Outside of Menu Buttons; track location for possible dragging
  if (mouseY > buttonHeight || menuIndex > 7)
  {
    prevMouseX = mouseX;
    prevMouseY = mouseY;
    return;
  }
  
  int id = menuIndex + 1;  // First Button, Second Button and so on
  if (id == 1 && menuID == id) return; 
  menuID = id;

  // Take action based on user selection
  switch (menuID) 
  {
    case 1:
      loadRandomBiologicalNetwork();
      break;

    case 2:
      loadCustomBiologicalNetwork();
      break;
      
    case 3:
      enrichPPIFromGeneSet();
      break;
    
    case 4:
      enrichPPIFromGeneExpression();
      break;
      
    case 5:
      determineHubProteins();
      break;
    
    case 6:
      zoomLevel = 1.0;
      panX = panY = 0;
      break;
    
    case 7:
      // Protein Structure Visualization. Not implemented yet
      break;
    
    case 8:
      // Shortest path. Handled when enter is pressed
      break;
    
    default:
      break;
  }  
}

/**
 * Key events handler
 */

void keyPressed() 
{
  if (key != CODED)
  {
    if (menuID == 8)
    {
      if (key == '\n' )
      { 
        savedUserText = typing; 
        typing        = ""; 
        determineShortestPath(); 
      } 
      else 
      { 
        typing = typing + key;  
        savedUserText = ""; 
      }
    }
    return;
  }
  
  if (keyCode == DOWN)  { panY += 50; return; }
  if (keyCode == UP)    { panY -= 50; return; }  
  if (keyCode == RIGHT) { panX += 50; return; }
  if (keyCode == LEFT)  { panX -= 50; return; }    
}
