#include <Eigen/Core>

class GravityFieldModel
{

public:
  // Constructor
  GravityFieldModel( double gravitationalParameter );

  // Function to calculate the potential
  double calculatePotential( Eigen::Vector3d relativePosition );

  // Function to calculate the potential
  Eigen::Vector3d calculatePotentialGradient( Eigen::Vector3d relativePosition );

private:
  // Property of object, determines strength of the field
  double gravitationalParameter_;
};

