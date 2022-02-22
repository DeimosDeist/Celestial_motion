import numpy as np

def potential(point1, point2):
    dv = 1/(np.absolute(point1-point2))

G = 10

def abs_value(vector):
    """
    Calculates the length of a vector

    Args:
        vector  (np.array)      Vector to get the length from the

    Returns:
        Length  (float)     Length of the Vector
    """
    return np.sqrt(np.sum(np.square(vector)))

class MassPoint:
    """

    """
    def __init__(self, position: np.array, mass: float):
        """
        Initializes a MassPoint with a position and a mass of

        Attr:
        position (np.array): Position of the MassPoints
        mass (float): mass of the MassPoint/Planet
        """
        self.__position = position
        self.__mass = mass


    @property
    def position(self):
        return self.__position

    @position.setter
    def position(self, value: np.array):
        self.__position = value

    @property
    def mass(self):
        return self.__mass


class MassSystem:
    """
    sets up a SpringMassSystem to be calculated and animated. Springs, Tethers and Masses need to be created beforehand!
    """

    def __init__(self, points, time_step=.001):
        """
        Initializes a SpringMassSystem that will be calculated
        Arguments:
        points (tuple): All Masses and Tethers of the MassSystem

        gravity (float): The gravity of the MassSystem. Pointing down = negative value. Default value is -9.8.
        time_step (float): the time steps to be integrated over. Default value is 0.0001, Make sure to not make it too small.
        points (tuple): All Masses and Tethers of the MassSystem
        springs (tuple): All Springs of the system
        k_matrix (np.array): all the spring stiffnesseseses in the system; the indices indicate which masspoints the springs connect

        Attr:
        self._points (tuple): All Masses and Tethers of the system
        self._gravity (float)= gravity that is acting on the MassSystem
        self._time_step (float)= time step that is used to integrate
        self._springs (tuple) = springs of the MassSystem


        """
        self._points = points
        self._time_step = time_step

    def __str__(self):
        """
        A nice text representation of the SpringMassSystem
        """
        return (
            f"We are looking at the planets {self.points}")

    @property
    def time_step(self):
        return self._time_step

    @property
    def points(self):
        return self._points

    def getforces(self):
        global G
        matrix = [['null' for j in range(len(self.points))] for i in range(len(self.points))]
        for i in range(len(self.points)):
            for j in range(len(self.points)):
                if i == j:
                    matrix[i][j] = np.array([0,0])
                    pass
                if not i == j:
                    matrix[i][j] = G*self.points[i].mass*self.points[j].mass*(self.points[i].position-self.points[j].position)/abs_value(self.points[i].position-self.points[j].position)**3
        return matrix

    def equation_of_motion(self, point, velocity, position):
        """
        Returns the equation of motion for the point


        Arguments: point : (MassPoint)         Masspoint whose acceleration is summed up over all other points (with
                                                the new_dv() function)
        velocity: (np.array)         Array containing the velocities
        position: (np.array)         Array containing the positions
        """

        dv = self.new_dv(point, position)
        dx = velocity[point, :]
        return dv, dx

    def new_dv(self, point, position):
        """
        Calculates the new dv with the spring matrix for point "point"

        Arguments:
            point : (MassPoint)
            position: (np.array)
        """
        for i in range(len(self.points)):
            if i != point:
                dv = dv + self.points[point].acceleration(position[point, :], position[i, :])
        return dv

    def change(self):
        new_vel = np.zeros((len(self.points), len(self.points[0].position)))
        new_pos = np.zeros((len(self.points), len(self.points[0].position)))
        k_1 = np.zeros((len(self.points), 2, len(self.points[0].position)))
        k_2 = np.zeros((len(self.points), 2, len(self.points[0].position)))
        k_3 = np.zeros((len(self.points), 2, len(self.points[0].position)))
        k_4 = np.zeros((len(self.points), 2, len(self.points[0].position)))

        for i in range(len(self.points)):
            new_vel[i, :] = self.points[i].velocity  # gets initial position into vel and pos (for the first step)
            new_pos[i, :] = self.points[i].position
        for i in range(len(self.points)):
            k_1[i, :, :] = self.equation_of_motion(i, new_vel, new_pos)
        for i in range(len(self.points)):
            k_2[i, :, :] = self.equation_of_motion(i, new_vel + (self.time_step / 2) * k_1[:, 0, :],
                                                   new_pos + (self.time_step / 2) * k_1[:, 1, :])
        for i in range(len(self.points)):
            k_3[i, :, :] = self.equation_of_motion(i, new_vel + (self.time_step / 2) * k_2[:, 0, :],
                                                   new_pos + (self.time_step / 2) * k_2[:, 1, :])
        for i in range(len(self.points)):
            k_4[i, :, :] = self.equation_of_motion(i, new_vel + self.time_step * k_3[:, 0, :],
                                                   new_pos + self.time_step * k_3[:, 1, :])
        for i in range(len(self.points)):
            if isinstance(self.points[i], Tether):
                new_vel[i, :] = 0
                new_pos[i, :] = self.points[i].position
            else:
                # new_vel[i,:] = self.points[i].velocity + dv*self.time_step
                new_vel[i, :] = self.points[i].velocity + self.time_step / 6 * (
                        k_1[i, 0, :] + 2 * k_2[i, 0, :] + 2 * k_3[i, 0, :] + k_4[i, 0, :])
                new_pos[i, :] = self.points[i].position + self.time_step / 6 * (
                        k_1[i, 1, :] + 2 * k_2[i, 1, :] + 2 * k_3[i, 1, :] + k_4[i, 1, :])
            self.points[i].position = new_pos[i, :]
            self.points[i].velocity = new_vel[i, :]

    def simulate(self, t):
        """
        Runs the simulation,

        Args:
            t:  (float)      the total time the simulation should run for.

        Returns:
            result, velocity:    (np.array,np.array)        Matrices with the result matrix and the velocity Matrix.
        """
        steps_taken = int(np.round(t / self.time_step))
        result = np.zeros((len(self.points[0].position), steps_taken + 1,
                           len(self.points)))  # contains [2 dimensions for space, the steps taken, and the points]

        velocity = np.zeros((len(self.points[1].position), steps_taken + 1, len(self.points)))
        for i in range(len(self.points)):
            result[:, 0, i] = self.points[i].position
            velocity[:, 0, i] = self.points[i].velocity
        for i in range(steps_taken):
            self.change()
            for j in range(len(self.points)):
                result[:, i + 1, j] = self.points[j].position
                velocity[:, i + 1, j] = self.points[j].velocity

        return result, velocity

    def energy_step(self, step, result_matrix, velocity_matrix):
        energy_mass_at_step = 0
        energy_spring_at_step = 0
        steps_taken = result_matrix.shape[1]

        for i in range(len(self.points)):
            v_squared = abs_value(velocity_matrix[:, step, i]) ** 2
            energy_mass_at_step = energy_mass_at_step + self.points[i].mass * -self.gravity[1] * (
            result_matrix[1, step, i]) + 0.5 * self.points[i].mass * v_squared

        for i in range(self.spring_matrix.shape[0]):
            for j in range(self.spring_matrix.shape[1]):
                spring_stretch = abs_value(result_matrix[:, step, i] - result_matrix[:, step, j])
                energy_spring_at_step = energy_spring_at_step + 0.25 * self.spring_matrix[i, j].stiffness * (
                        self.spring_matrix[i, j].rest_length - spring_stretch) ** 2

        return energy_spring_at_step + energy_mass_at_step

    def energy(self, result_matrix, velocity_matrix):
        steps_taken = result_matrix.shape[1]
        energy = np.zeros(steps_taken)
        for i in range(steps_taken):
            energy[i] = self.energy_step(i, result_matrix, velocity_matrix)
        return energy

    def kin_energy_step(self, step, velocity_matrix):
        kin_energy_at_step = 0

        for i in range(len(self.points)):
            v_squared = abs_value(velocity_matrix[:, step, i]) ** 2
            kin_energy_at_step = kin_energy_at_step + 1 / 2 * self.points[i].mass * v_squared
        return kin_energy_at_step

    def kin_energy(self, result_matrix, velocity_matrix):
        steps_taken = result_matrix.shape[1]
        kin_energy = np.zeros(steps_taken)
        for i in range(steps_taken):
            kin_energy[i] = self.kin_energy_step(i, velocity_matrix)
        return kin_energy