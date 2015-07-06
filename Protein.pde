/**
 * Basic Protein class. Proteins are nodes of the Protein-Protein 
 * Interaction (PPI) network. The class uses Processing 2.0 class 
 * PVector for (x, y, z) location as well as current and previous
 * node velocities. The protein velocity as well as other parameters
 * are used in force layout simulation for elegant display of
 * biological network network.
 *
 * The instance variable 'links' keeps track of all interactions.
 *
 * @author Pranav Srinivas
 * @version 12 May 2015
 */

class Protein
{
  private String name;                  // Protein Name (Gene Symbol)
  private String id;                    // Unique protein ID
  private float xMin, yMin, zMin;       // Boundary positions
  private float xMax, yMax, zMax;       // Boundary positions
  private float velocityLimit;          // Maximum velocity for force simulation
  private float repulsionStrength;      // -ve: Repulsion, +ve Attraction  
  private PVector location;             // Position of this protein
  private PVector currentVelocity;      // Current protein velocity
  private PVector previousVelocity;     // Last protein velocity
  private float   radius;               // Influence radius of protein
  private float   dampingFactor;        // Force simulation damping factor 
  private float   forceRampFactor;      // Force simulation ramp factor
  private float   diameter;             // Diameter for drawing
  private boolean isEnriched;           // A boolean flag for enrichment
  private boolean isHub;                // Is Hub Protein ?
  private int     linksCount;           // Enriched link count for enrichment
  private boolean freeze;               // Freeze protein movement under certain conditions
  private boolean isVisited;            // flag for use in shortest path
  private ArrayList<Interaction> links; // All links of this protein (used in path queries)
  
  /**
   * Initialize data members. Helper method for constructor
   */
  void init()
  {
    name     = "";
    id       = "";
    radius   = 250;
    diameter = 0;

    xMin = yMin = zMin = -Float.MAX_VALUE;
    xMax = yMax = zMax =  Float.MAX_VALUE;
    velocityLimit      = 10;
    dampingFactor      = 0.5;
    repulsionStrength  = -1;
    forceRampFactor    = 1.0;

    location         = new PVector(0, 0, 0);
    currentVelocity  = new PVector();
    previousVelocity = new PVector();
    isEnriched       = false;
    isHub            = false;
    linksCount       = 0;
    freeze           = false;
    isVisited        = false;
    links            = new ArrayList<Interaction>();
  }
  
  /** 
   * Assign a random location
   */
  void initLocation()
  {
    location.set( width/2+random(-0.4*width, 0.4*width), height/2+random(-0.4*height, 0.4*height) );
    freeze = false;
  } 

  //-------------------------------------------------------------
  // Protein Constructors
  //-------------------------------------------------------------

  Protein()
  {
    super();
    init();
  }

  Protein(float X, float Y)
  {
    super();
    init();
    location.set( X, Y );
  }

  Protein(PVector vector)
  {
    super();
    init();
    location.set( vector );
  }

  Protein(float X, float Y, float Z)
  {
    super();
    init();
    location.set( X, Y, Z );
  }
  
  Protein(String pname)
  {
    super();
    init();
    initLocation();
    setName(pname);
  }

  //-------------------------------------------------------------
  // Accessor Methods (Getters)
  //-------------------------------------------------------------

  // Get X location of Protein
  float getX() 
  {
    return location.x;
  }

  // Get Y location of Protein
  float getY()
  {
    return location.y;
  }

  // Get Z location of Protein
  float getZ()
  {
    return location.z;
  }

  String getName()
  {
    return name;
  }
  
  String getId()
  {
    return id;
  }

  // Get force ramp factor
  float getForceRampFactor()
  {
    return forceRampFactor;
  }

  // Get damping Factor
  float getDampingFactor()
  {
    return dampingFactor;
  }

  // Get radius of influence
  float getRadius()
  {
    return radius;
  }

  float getDiameter()
  {
    return diameter;
  }

  float getVelocityLimit()
  {
    return velocityLimit;
  }

  float getRepulsionStrength()
  {
    return repulsionStrength;
  }

  PVector getLocation()
  {
    return location.get();
  }
  
  boolean getIsEnriched()
  {
    return isEnriched;
  }
  
  boolean getIsHub()
  {
    return isHub;
  }

  int getLinksCount()
  {
    return linksCount;
  }
  
  boolean getIsVisited()
  {
    return isVisited;
  }
  
  ArrayList<Interaction> getLinks()
  {
    return links;
  }

  //-------------------------------------------------------------
  // Mutator Methods (Setters)
  //-------------------------------------------------------------

  // Set X location of Protein
  void setX(float x) 
  {
    location.x = x;
  }

  // Set Y location of Protein
  void setY(float y)
  {
    location.y = y;
  }

  // Set Z location of Protein
  void setZ(float z)
  {
    location.z = z;
  }

  // Add to X location of Protein
  void addX(float x) 
  {
    location.x += x;
  }

  // Add to Y location of Protein
  void addY(float y)
  {
    location.y += y;
  }

  // Add to Z location of Protein
  void addZ(float z)
  {
    location.z += z;
  }

  // Set force ramp factor
  void setForceRampFactor(float ramp)
  {
    forceRampFactor = ramp;
  }

  // Set damping Factor
  void setDampingFactor(float damp)
  {
    dampingFactor = damp;
  }

  void setName(String pname)
  {
    name = pname;
  }

  void setId(String pid)
  {
    id = pid;
  }

  void setRadius(float r)
  {
    radius = r;
  }

  void setDiameter(float d)
  {
    diameter = d;
  }

  /**
   * Set boundary locations
   * @param MinX, MinY, MinZ, MaxX, MaxY, MaxZ (boundary positions)
   */
  void setPositionLimit(float MinX, float MinY, float MinZ, float MaxX, float MaxY, float MaxZ)
  {
    xMin = MinX; xMax = MaxX;
    yMin = MinY; yMax = MaxY;
    zMin = MinZ; zMax = MaxZ;
  }

  /** 
   * Set boundary locations
   * @param MinX, MinY, MaxX, MaxY (boundary positions)
   */
  void setPositionLimit(float MinX, float MinY, float MaxX, float MaxY)
  {
    xMin = MinX; xMax = MaxX;
    yMin = MinY; yMax = MaxY;
  }

  void setVelocityLimit(float vl)
  {
    velocityLimit = vl;
  }

  void setRepulsionStrength(float strength)
  {
    repulsionStrength = strength;
  }

  /**
   * Add delta velocity vector to current velocity vector of the protein
   */
  void addToCurrentVelocity(PVector delta)
  {
    currentVelocity.add(delta);
  }

  /** 
   * Set or clear isEnriched boolean flag
   */
  void setIsEnriched(boolean value)
  {
    isEnriched = value;
  }
  
  /** 
   * Mark/Unmark as hub
   */
  void setIsHub(boolean value)
  {
    isHub = value;
    if (value == true) setRadius(140);
    else               setRadius(100);
  }
  
  void setIsVisited(boolean value)
  {
    isVisited = value;
  }

  void incrementLinksCount()
  {
    ++linksCount;
  }

  void clearLinksCount()
  {
    linksCount = 0;
  }
  
  void addLink(Interaction link)
  {
    links.add(link);
  }
  
  void clearLinks()
  {
    links.clear();
  }
  
  /** 
   * Reset location and velocity
   */
  void reset()
  {
    location.set(width/2+random(-200, 200), height/2+random(-200, 200));
    currentVelocity.set(0, 0, 0);
    previousVelocity.set(0, 0, 0); 
    isEnriched = false;
    freeze = false;
    setIsHub(false);
  }

  /**
   * Protein Repulsion and Attraction 
   * Repulse neighbor protein 'nbr' based on distance and radius
   * @parm nbr neighboring protein node
   */
  void repulse(Protein nbr)
  {
    float distance = PVector.dist(this.getLocation(), nbr.getLocation());

    // If 'nbr' is too close then don't exert any force
    if (distance <= 0) return;

    // If 'nbr' is too far away then don't exert any force
    if (distance >= radius) return;

    // Force formula adapted from http://cs.brown.edu/~rt/gdhandbook/chapters/force-directed.pdf

    float p = pow(distance/radius, 1/forceRampFactor);
    float force = repulsionStrength * p * 9 * (1 / (p + 1) + ((p - 3) / 4)) / distance;

    // Calculate Delta force vector as difference vector times force
    PVector deltaForce = PVector.sub(this.getLocation(), nbr.getLocation());
    deltaForce.mult(force);

    // Change current velocity of neighbor protein 'nbr'
    nbr.currentVelocity.add(deltaForce.x, deltaForce.y, deltaForce.z);
  }

  /**
   * Repulse ArrayList of neighbors
   * @param neighbors ArrayList of neighboring proteins 
   */
  void repulse(ArrayList<Protein> neighbors)
  {
    for (Protein nbr : neighbors)
    {
      if (nbr != this) repulse(nbr); // Repulse everyone except self!
    }
  }

  //-------------------------------------------------------------
  // Rotation Methods 
  //-------------------------------------------------------------

  /** 
   * Rotate radian angle around X axis
   * @param radian angle of rotation
   */
  void rotateX(float radian)
  {
    float nextY = getY() * cos(radian) - getZ() * sin(radian);
    float nextZ = getY() * sin(radian) + getZ() * cos(radian);
    setY( nextY );
    setZ( nextZ );
  }

  /** 
   * Rotate radian angle around Y axis
   * @param radian angle of rotation
   */
  void rotateY(float radian)
  {
    float nextX = getX() * cos(-radian) - getZ() * sin(-radian);
    float nextZ = getX() * sin(-radian) + getZ() * cos(-radian);
    setX( nextX );
    setZ( nextZ );
  }

  /**
   * Rotate radian angle around Z axis
   * @param radian angle of rotation
   */
  void rotateZ(float radian)
  {
    float nextX = getX() * cos(radian) - getY() * sin(radian);
    float nextY = getX() * sin(radian) + getY() * cos(radian);
    setX( nextX );
    setY( nextY );
  }

  /**
   * Update protein velocity and location
   * @param updateX update x component of position
   * @param updateY update y component of position
   * @param updateZ update z component of position
   */
  void updateVelocityAndPosition(boolean updateX, boolean updateY, boolean updateZ)
  {
    // Copy current velocity into previous velocity before the update
    previousVelocity = currentVelocity.get();

    // Update location based on current velocity
    if (updateX) addX( currentVelocity.x );
    if (updateY) addY( currentVelocity.y );
    if (updateZ) addZ( currentVelocity.z );

    // Check for out of range locations
    handleBoundaryConditions();

    // Update current velocity for damping
    currentVelocity.mult(1 - dampingFactor);

    // Limit the current velocity magnitude to specified value
    currentVelocity.limit( velocityLimit );
  }

  /** 
   * Update velocity and protein position
   */
  void updateVelocityAndPosition()
  {
    if (freeze) return;  // Frozen proteins don't move!
    updateVelocityAndPosition(true, true, true);
  }

  /** 
   * Adjust position and reverse currentVelocity around the edges
   */
  void handleBoundaryConditions()
  {
    float x = getX();
    float y = getY();
    float z = getZ();

    if (x < xMin)
    {
      x = xMin + random(400, 500);
      currentVelocity.x = -currentVelocity.x;
      freeze = true;
    }
    if (x > xMax)
    {
      x = xMax - random(400, 500);
      currentVelocity.x = -currentVelocity.x;
      freeze = true;
    }

    if (y < yMin) 
    {
      y = yMin + random(300, 400);
      currentVelocity.y = -currentVelocity.y;
      freeze = true;
    }
    if (y > yMax) 
    {
      y = yMax - random(300, 400);
      currentVelocity.y = -currentVelocity.y;
      freeze = true;
    }

    if (z < zMin) 
    {
      z = zMin + random(300, 400);
      currentVelocity.z = -currentVelocity.z;
      freeze = true;
    }
    if (z > zMax) 
    {
      z = zMax - random(300, 400);
      currentVelocity.z = -currentVelocity.z;
      freeze = true;
    }

    // Update location for revised values of x, y and z
    location.set(x, y, z);
  }

}
