import numpy as np

class MassPoint:
    """
    Masspoint is a point with a mass that is fixed to a spring.
    """

    def __init__(self, mass, position, velocity):
        """
        Construct a MassPoint

        Attr:
        self._mass = Mass of the MassPoint
        self._position(np.array) = Position of the MassPoint
        self._velocity(np.array) = Velocity of the MassPoint
        """
        self._mass = mass
        self._position = position
        self._velocity = velocity

    def __str__(self):
        return str(self._mass + " at " + self._position + " with speed " + self._velocity)

    # setters and getters for MassPoint

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, value):
        self._mass = value

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, value):
        self._position = value

    @property
    def velocity(self):
        return self._velocity

    @velocity.setter
    def velocity(self, value):
        self._velocity = value

    def acceleration(self, spring, position0, position1):
        """
        gives acceleration for 1 point at 1 spring
        Args:
            spring      (Spring)    Spring for which to calculate the effect on the mass point
            position0   (np.arrays)    Position of the spring
            position1   (np.arrays)    Position of the spring

        Returns:
        dv (np.array): acceleration
        """
        if abs_value(position0 - position1) == 0:
            dv = 0
        else:
            dv = spring.stiffness * (abs_value(position0 - position1) - spring.rest_length) / self.mass * (
                    position1 - position0) / (abs_value(position0 - position1))
        return dv
