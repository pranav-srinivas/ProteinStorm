/**
 * Basic Interaction class. Interactions are links between two Proteins
 * In addition to the source and sink proteins the class also has
 * adjustable parameters (rigidity, length, dampingFactor) for
 * force layout simulation.
 *
 * @author Pranav Srinivas
 * @version 14 May 2015
 */

class Interaction
{
  private Protein sourceProtein;      // Source protein
  private Protein sinkProtein;        // Sink protein
  private float dampingFactor;        // Force damping factor
  private float length;               // length
  private float rigidity;             // Interaction rigidity (strength)
  private boolean isVisited;          // flag used in shortest path query
  private boolean highlight;          // flag used in shortest path visualization

  //-------------------------------------------------------------
  // Interaction Constructors
  //-------------------------------------------------------------

  Interaction(Protein startProtein, Protein endProtein)
  {
    sourceProtein = startProtein;
    sinkProtein   = endProtein;
    dampingFactor = 0.85;
    length        = 100;
    rigidity      = 0.7;
    highlight     = false;
    isVisited     = false;
  }

  Interaction(Protein startProtein, Protein endProtein, float damping, float Length, float firmness)
  {
    sourceProtein = startProtein;
    sinkProtein   = endProtein;
    dampingFactor = damping;
    length        = Length;
    rigidity      = firmness;
  }

  //-------------------------------------------------------------
  // Interaction Accessors (Getters)
  //-------------------------------------------------------------

  Protein getSourceProtein()
  {
    return sourceProtein;
  }

  Protein getSinkProtein()
  {
    return sinkProtein;
  }

  float getDampingFactor()
  {
    return dampingFactor;
  }

  float getLength()
  {
    return length;
  }

  float getRigidity()
  {
    return rigidity;
  }
  
  boolean getHighlight()
  {
    return highlight;
  }
    
  boolean getIsVisited()
  {
    return isVisited;
  }
  
  //-------------------------------------------------------------
  // Interaction Mutators (Setters)
  //-------------------------------------------------------------

  void setSourceProtein(Protein protein)
  {
    sourceProtein = protein;
  }

  void setSinkProtein(Protein protein)
  {
    sinkProtein = protein;
  }

  void setDampingFactor(float df)
  {
    dampingFactor = df;
  }

  void setLength(float l)
  {
    length = l;
  }

  void setRigidity(float r)
  {
    rigidity = r;
  }

  void setHighlight(boolean value)
  {
    highlight = value;
  }
  
  void setIsVisited(boolean value)
  {
    isVisited = value;
  }
  
  /**
   * Apply forces to source and sink proteins
   */
  void applyForces()
  {
    // Get delta vector of location vectors. Normalize and scale by length
    PVector delta = PVector.sub(sinkProtein.getLocation(), sourceProtein.getLocation());
    delta.normalize();
    delta.mult(length);

    // Derive target position vector by adding the delta vector to from protein's location
    PVector target = PVector.add(sourceProtein.getLocation(), delta);

    // Calculate delta force vectors for to and from proteins.
    // Newton's 3rd law: sourceForce = -sinkForce
    PVector sinkForce = PVector.sub(target, sinkProtein.getLocation());
    sinkForce.mult(0.5);
    sinkForce.mult(rigidity);
    sinkForce.mult(1 - dampingFactor);
    PVector sourceForce = PVector.mult(sinkForce, -1);

    sinkProtein.addToCurrentVelocity( sinkForce );
    sourceProtein.addToCurrentVelocity( sourceForce );
  }

}
