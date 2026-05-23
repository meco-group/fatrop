import casadi as ca
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, interp1d

class CarTrajectoryOptimizer:
    def __init__(self):
        # Vehicle parameters (typical car values)
        self.L = 2.7  # wheelbase [m]
        self.max_speed = 20.0  # maximum speed [m/s]
        self.max_acceleration = 3.0  # maximum acceleration [m/s^2]
        self.max_steering = np.pi/4  # maximum steering angle [rad]
        
        # Optimization parameters
        self.N = 50  # number of control intervals (K-1 in fatrop notation)
        self.K = self.N + 1  # number of time steps (K in fatrop notation)
        self.T = 5.0  # time horizon [s]
        self.dt = self.T / self.N  # time step
        
        # Reference tracking weights
        self.Q_pos = 10.0  # position tracking weight
        self.Q_heading = 10.  # heading tracking weight
        self.Q_velocity = .1  # velocity tracking weight
        self.R_acc = 1.  # acceleration penalty weight
        self.R_steer = 0.0  # steering penalty weight
        
        # Define dimensions following fatrop structure
        self.nx = [4 for _ in range(self.K)]  # state dimensions
        self.nu = [2 for _ in range(self.N)] + [0]  # control dimensions
        
        # Reference trajectory storage
        self.reference_trajectory = None
    
    def generate_spline_reference(self, waypoints, desired_speed=5.0):
        """
        Generate a spline reference trajectory from waypoints
        
        Args:
            waypoints: List of (x, y) waypoints
            desired_speed: Desired speed along the trajectory [m/s]
        
        Returns:
            Dictionary containing spline interpolators and reference data
        """
        waypoints = np.array(waypoints)
        
        # Calculate cumulative distance along waypoints
        distances = np.zeros(len(waypoints))
        for i in range(1, len(waypoints)):
            distances[i] = distances[i-1] + np.linalg.norm(waypoints[i] - waypoints[i-1])
        
        # Create splines for x and y coordinates
        spline_x = CubicSpline(distances, waypoints[:, 0])
        spline_y = CubicSpline(distances, waypoints[:, 1])
        
        # Calculate total path length
        total_length = distances[-1]
        
        # Generate reference trajectory at optimization time points
        time_points = np.linspace(0, self.T, self.K)
        
        # Map time to distance along path (assuming constant speed)
        path_distances = (desired_speed * time_points) % total_length
        
        # Evaluate splines at these distances
        x_ref = spline_x(path_distances)
        y_ref = spline_y(path_distances)
        
        # Calculate reference heading (tangent to path)
        dx_ds = spline_x.derivative()(path_distances)
        dy_ds = spline_y.derivative()(path_distances)
        theta_ref = np.arctan2(dy_ds, dx_ds)
        
        # Set reference velocity
        v_ref = np.full(self.K, desired_speed)
        
        # Store reference trajectory
        self.reference_trajectory = {
            'x': x_ref,
            'y': y_ref,
            'theta': theta_ref,
            'v': v_ref,
            'spline_x': spline_x,
            'spline_y': spline_y,
            'total_length': total_length,
            'desired_speed': desired_speed
        }
        
        return self.reference_trajectory
    
    def evaluate_reference_at_time(self, t):
        """
        Evaluate reference trajectory at a specific time
        
        Args:
            t: Time [s]
        
        Returns:
            Reference state [x, y, theta, v]
        """
        if self.reference_trajectory is None:
            raise ValueError("No reference trajectory set. Call generate_spline_reference first.")
        
        # Map time to distance along path
        distance = (self.reference_trajectory['desired_speed'] * t) % self.reference_trajectory['total_length']
        
        # Evaluate splines
        x_ref = self.reference_trajectory['spline_x'](distance)
        y_ref = self.reference_trajectory['spline_y'](distance)
        
        # Calculate heading
        dx_ds = self.reference_trajectory['spline_x'].derivative()(distance)
        dy_ds = self.reference_trajectory['spline_y'].derivative()(distance)
        theta_ref = np.arctan2(dy_ds, dx_ds)
        
        v_ref = self.reference_trajectory['desired_speed']
        
        return np.array([x_ref, y_ref, theta_ref, v_ref])
        
    def discrete_dynamics(self, uk, xk, k):
        """
        Discrete dynamics using RK4 integration
        Following fatrop structure: returns x_{k+1}
        """
        # RK4 integration of bicycle model
        def bicycle_model_continuous(x, u):
            x_pos, y_pos, theta, v = x[0], x[1], x[2], x[3]
            a, delta = u[0], u[1]
            
            x_dot = v * ca.cos(theta)
            y_dot = v * ca.sin(theta)
            theta_dot = v * ca.tan(delta) / self.L
            v_dot = a
            
            return ca.vertcat(x_dot, y_dot, theta_dot, v_dot)
        
        k1 = bicycle_model_continuous(xk, uk)
        k2 = bicycle_model_continuous(xk + self.dt/2 * k1, uk)
        k3 = bicycle_model_continuous(xk + self.dt/2 * k2, uk)
        k4 = bicycle_model_continuous(xk + self.dt * k3, uk)
        
        return xk + self.dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    
    def cost(self, uk, xk, k):
        """
        Stage cost function for reference tracking
        """
        cost_val = 0
        
        # Reference tracking costs (if reference trajectory is available)
        if self.reference_trajectory is not None:
            # Get reference state at this time step
            x_ref = self.reference_trajectory['x'][k]
            y_ref = self.reference_trajectory['y'][k]
            theta_ref = self.reference_trajectory['theta'][k]
            v_ref = self.reference_trajectory['v'][k]
            
            # Position tracking error
            pos_error = (xk[0] - x_ref)**2 + (xk[1] - y_ref)**2
            cost_val += self.Q_pos * pos_error
            
            # Heading tracking error (handle angle wrapping)
            heading_error = ca.sin(xk[2] - theta_ref)**2 + (1 - ca.cos(xk[2] - theta_ref))**2
            cost_val += self.Q_heading * heading_error
            
            # Velocity tracking error
            velocity_error = (xk[3] - v_ref)**2
            cost_val += self.Q_velocity * velocity_error
        
        # Control effort penalties
        if k < self.N:
            cost_val += self.R_acc * uk[0]**2  # acceleration penalty
            cost_val += self.R_steer * uk[1]**2  # steering penalty
        
        return cost_val
    
    def path_constraints(self, uk, xk, k, start_state, goal_state=None, obstacles=None):
        """
        Path constraints following fatrop structure
        For reference tracking, goal_state can be None (no terminal constraints)
        """
        cc = []
        
        # Equality constraints
        if k == 0:
            # Initial condition
            cc.append(xk[0] - start_state[0] == 0.)
            cc.append(xk[1] - start_state[1] == 0.)
            cc.append(xk[2] - start_state[2] == 0.)
            cc.append(xk[3] - start_state[3] == 0.)
        elif k == self.K - 1 and goal_state is not None and self.reference_trajectory is None:
            # Terminal condition (only for point-to-point, not reference tracking)
            cc.append(xk[0] - goal_state[0] == 0.)
            cc.append(xk[1] - goal_state[1] == 0.)
            cc.append(xk[2] - goal_state[2] == 0.)
            cc.append(xk[3] - goal_state[3] == 0.)
        
        # Inequality constraints
        # State bounds
        if k > 0:
            cc.append(-self.max_speed / 4 <= (xk[3] <= self.max_speed))  # velocity bounds
        
        # Control bounds (only if controls exist)
        if k < self.N:
            cc.append(-self.max_acceleration <= (uk[0] <= self.max_acceleration))  # acceleration bounds
            cc.append(-self.max_steering <= (uk[1] <= self.max_steering))  # steering bounds
        
        # Obstacle avoidance constraints
        if obstacles is not None:
            for obs in obstacles:
                obs_x, obs_y, obs_radius = obs
                dist = (xk[0] - obs_x)**2 + (xk[1] - obs_y)**2
                cc.append(obs_radius**2 <= dist)  # minimum clearance
        
        return cc
    
    def setup_optimization_problem(self, start_state, goal_state, obstacles=None):
        """
        Set up the trajectory optimization problem following fatrop structure
        """
        # Create optimization problem
        opti = ca.Opti()
        
        # Decision variables - following fatrop structure
        x = []
        u = []
        ng = []
        for k in range(self.K):
            x.append(opti.variable(self.nx[k]))
            u.append(opti.variable(self.nu[k]))
        
        # Add constraints - following fatrop order
        for k in range(self.K):
            # Dynamics constraints
            if k < self.K - 1:
                opti.subject_to(x[k+1] == self.discrete_dynamics(u[k], x[k], k))
            
            # Path constraints
            path_constr = self.path_constraints(u[k], x[k], k, start_state, goal_state, obstacles)
            for constr in path_constr:
                opti.subject_to(constr)
            ng.append(sum([ci.nnz() for ci in path_constr]))  # number of constraints for this step
        
        # Set the objective - following fatrop structure
        J = 0
        for k in range(self.K):
            J += self.cost(u[k], x[k], k)
        
        opti.minimize(J)
        
        # Initial guess - use reference trajectory if available, otherwise linear interpolation
        for k in range(self.K):
            if self.reference_trajectory is not None:
                # Use reference trajectory as initial guess
                opti.set_initial(x[k], ca.vertcat(
                    self.reference_trajectory['x'][k],
                    self.reference_trajectory['y'][k],
                    self.reference_trajectory['theta'][k],
                    self.reference_trajectory['v'][k]
                ))
            elif goal_state is not None:
                # Linear interpolation between start and goal
                for i in range(4):
                    opti.set_initial(x[k][i], start_state[i] + (goal_state[i] - start_state[i]) * k / self.N)
            else:
                # Use start state for all time steps if no goal specified
                for i in range(4):
                    opti.set_initial(x[k][i], start_state[i])
            
            # Zero initial guess for controls
            if k < self.N:
                opti.set_initial(u[k], ca.vertcat(0, 0))
        
        # Solver options - try fatrop first, fallback to ipopt
        opti.solver('fatrop', {
            'structure_detection': 'manual', 
            'nx': self.nx, 
            'nu': self.nu, 
            'ng': ng, 
            'N': self.N, 
            "expand": True, 
            "fatrop.tol": 1e-6,
            "jit": False
        })
        # ipopt with lbfgs and a large memory size
        # opti.solver("ipopt", {"ipopt.hessian_approximation": "limited-memory", "ipopt.tol": 1e-6})
        print("Using FATROP solver")
        
        return opti, x, u
    
    def solve_trajectory(self, start_state, goal_state, obstacles=None, warm_start_data=None):
        """
        Solve the trajectory optimization problem
        """
        print("Setting up optimization problem using Opti stack...")
        opti, x, u = self.setup_optimization_problem(start_state, goal_state, obstacles)
        
        # Apply warm start if provided
        if warm_start_data is not None:
            print("Applying warm start from previous solution...")
            prev_opti, prev_x, prev_u, prev_sol = warm_start_data
            
            # Set initial guess using previous solution
            for k in range(self.K):
                opti.set_initial(x[k], prev_sol.value(prev_x[k]))
                if k < self.N:
                    opti.set_initial(u[k], prev_sol.value(prev_u[k]))
        
        print("Solving optimization problem...")
        try:
            # Solve the problem
            sol = opti.solve()
            
            # Extract solution - convert from list format to matrix format for compatibility
            X_opt = np.zeros((4, self.K))
            U_opt = np.zeros((2, self.N))
            
            for k in range(self.K):
                X_opt[:, k] = sol.value(x[k])
                if k < self.N:
                    U_opt[:, k] = sol.value(u[k])
            
            cost = sol.value(opti.f)
            
            return X_opt, U_opt, cost, True, (opti, x, u, sol)
            
        except Exception as e:
            print(f"Optimization failed: {e}")
            # Try to get debug solution
            try:
                X_opt = np.zeros((4, self.K))
                U_opt = np.zeros((2, self.N))
                
                for k in range(self.K):
                    X_opt[:, k] = opti.debug.value(x[k])
                    if k < self.N:
                        U_opt[:, k] = opti.debug.value(u[k])
                
                cost = opti.debug.value(opti.f)
                return X_opt, U_opt, cost, False, None
            except:
                return None, None, None, False, None
    
    def plot_results(self, X_opt, U_opt, obstacles=None, success=True):
        """
        Plot the optimization results
        """
        time = np.linspace(0, self.T, self.N + 1)
        time_control = np.linspace(0, self.T - self.dt, self.N)
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        
        # Set title based on success
        title_suffix = " (Converged)" if success else " (Debug - Not Converged)"
        
        # Trajectory plot
        color = 'b-' if success else 'r--'
        axes[0, 0].plot(X_opt[0, :], X_opt[1, :], color, linewidth=2, label='Actual Trajectory')
        axes[0, 0].plot(X_opt[0, 0], X_opt[1, 0], 'go', markersize=10, label='Start')
        axes[0, 0].plot(X_opt[0, -1], X_opt[1, -1], 'ro', markersize=10, label='End')
        
        # Plot reference trajectory if available
        if self.reference_trajectory is not None:
            axes[0, 0].plot(self.reference_trajectory['x'], self.reference_trajectory['y'], 
                          'k--', linewidth=2, alpha=0.7, label='Reference Trajectory')
        
        # Plot obstacles if provided
        if obstacles is not None:
            for obs in obstacles:
                circle = plt.Circle((obs[0], obs[1]), obs[2], color='red', alpha=0.3)
                axes[0, 0].add_patch(circle)
                # Safety margin
                circle_safe = plt.Circle((obs[0], obs[1]), obs[2] + 0.5, 
                                       color='orange', alpha=0.2, linestyle='--', fill=False)
                axes[0, 0].add_patch(circle_safe)
        
        axes[0, 0].set_xlabel('X [m]')
        axes[0, 0].set_ylabel('Y [m]')
        axes[0, 0].set_title('Vehicle Trajectory' + title_suffix)
        axes[0, 0].legend()
        axes[0, 0].grid(True)
        axes[0, 0].axis('equal')
        
        # Velocity profile
        axes[0, 1].plot(time, X_opt[3, :], color, linewidth=2, label='Actual Velocity')
        if self.reference_trajectory is not None:
            axes[0, 1].plot(time, self.reference_trajectory['v'], 'k--', linewidth=2, alpha=0.7, label='Reference Velocity')
        axes[0, 1].axhline(y=self.max_speed, color='k', linestyle='--', alpha=0.5, label='Max speed')
        axes[0, 1].set_xlabel('Time [s]')
        axes[0, 1].set_ylabel('Velocity [m/s]')
        axes[0, 1].set_title('Velocity Profile' + title_suffix)
        axes[0, 1].legend()
        axes[0, 1].grid(True)
        
        # Heading angle
        axes[0, 2].plot(time, X_opt[2, :] * 180/np.pi, color, linewidth=2, label='Actual Heading')
        if self.reference_trajectory is not None:
            axes[0, 2].plot(time, self.reference_trajectory['theta'] * 180/np.pi, 'k--', linewidth=2, alpha=0.7, label='Reference Heading')
        axes[0, 2].set_xlabel('Time [s]')
        axes[0, 2].set_ylabel('Heading [deg]')
        axes[0, 2].set_title('Heading Angle' + title_suffix)
        axes[0, 2].legend()
        axes[0, 2].grid(True)
        
        # Acceleration
        control_color = 'r-' if success else 'r--'
        axes[1, 0].plot(time_control, U_opt[0, :], control_color, linewidth=2)
        axes[1, 0].axhline(y=self.max_acceleration, color='k', linestyle='--', alpha=0.5)
        axes[1, 0].axhline(y=-self.max_acceleration, color='k', linestyle='--', alpha=0.5)
        axes[1, 0].set_xlabel('Time [s]')
        axes[1, 0].set_ylabel('Acceleration [m/s²]')
        axes[1, 0].set_title('Acceleration Input' + title_suffix)
        axes[1, 0].grid(True)
        
        # Steering angle
        axes[1, 1].plot(time_control, U_opt[1, :] * 180/np.pi, control_color, linewidth=2)
        axes[1, 1].axhline(y=self.max_steering * 180/np.pi, color='k', linestyle='--', alpha=0.5)
        axes[1, 1].axhline(y=-self.max_steering * 180/np.pi, color='k', linestyle='--', alpha=0.5)
        axes[1, 1].set_xlabel('Time [s]')
        axes[1, 1].set_ylabel('Steering Angle [deg]')
        axes[1, 1].set_title('Steering Input' + title_suffix)
        axes[1, 1].grid(True)
        
        # Tracking errors (if reference trajectory exists)
        if self.reference_trajectory is not None:
            pos_error = np.sqrt((X_opt[0, :] - self.reference_trajectory['x'])**2 + 
                              (X_opt[1, :] - self.reference_trajectory['y'])**2)
            axes[1, 2].plot(time, pos_error, 'g-', linewidth=2, label='Position Error')
            
            heading_error = np.abs(X_opt[2, :] - self.reference_trajectory['theta']) * 180/np.pi
            # Handle angle wrapping
            heading_error = np.minimum(heading_error, 360 - heading_error)
            axes[1, 2].plot(time, heading_error, 'r-', linewidth=2, label='Heading Error')
            
            axes[1, 2].set_xlabel('Time [s]')
            axes[1, 2].set_ylabel('Error')
            axes[1, 2].set_title('Tracking Errors' + title_suffix)
            axes[1, 2].legend()
            axes[1, 2].grid(True)
        else:
            # Phase portrait (x-y with velocity color coding) if no reference
            scatter = axes[1, 2].scatter(X_opt[0, :], X_opt[1, :], c=X_opt[3, :], 
                                       cmap='viridis', s=50)
            axes[1, 2].set_xlabel('X [m]')
            axes[1, 2].set_ylabel('Y [m]')
            axes[1, 2].set_title('Trajectory with Velocity' + title_suffix)
            plt.colorbar(scatter, ax=axes[1, 2], label='Velocity [m/s]')
            axes[1, 2].grid(True)
        
        plt.tight_layout()
        plt.show()
    
    def solve_reference_tracking(self, start_state, obstacles=None, warm_start_data=None):
        """
        Solve the reference tracking optimization problem
        """
        if self.reference_trajectory is None:
            raise ValueError("No reference trajectory set. Call generate_spline_reference first.")
        
        print("Setting up reference tracking optimization problem...")
        opti, x, u = self.setup_optimization_problem(start_state, goal_state=None, obstacles=obstacles)
        
        # Apply warm start if provided
        if warm_start_data is not None:
            print("Applying warm start from previous solution...")
            prev_opti, prev_x, prev_u, prev_sol = warm_start_data
            
            # Set initial guess using previous solution
            for k in range(self.K):
                opti.set_initial(x[k], prev_sol.value(prev_x[k]))
                if k < self.N:
                    opti.set_initial(u[k], prev_sol.value(prev_u[k]))
        else:
            # Use reference trajectory as initial guess
            for k in range(self.K):
                opti.set_initial(x[k], ca.vertcat(
                    self.reference_trajectory['x'][k],
                    self.reference_trajectory['y'][k],
                    self.reference_trajectory['theta'][k],
                    self.reference_trajectory['v'][k]
                ))
                if k < self.N:
                    opti.set_initial(u[k], ca.vertcat(0, 0))
        
        print("Solving reference tracking optimization problem...")
        try:
            # Solve the problem
            sol = opti.solve()
            
            # Extract solution
            X_opt = np.zeros((4, self.K))
            U_opt = np.zeros((2, self.N))
            
            for k in range(self.K):
                X_opt[:, k] = sol.value(x[k])
                if k < self.N:
                    U_opt[:, k] = sol.value(u[k])
            
            cost = sol.value(opti.f)
            
            return X_opt, U_opt, cost, True, (opti, x, u, sol)
            
        except Exception as e:
            print(f"Optimization failed: {e}")
            # Try to get debug solution
            try:
                X_opt = np.zeros((4, self.K))
                U_opt = np.zeros((2, self.N))
                
                for k in range(self.K):
                    X_opt[:, k] = opti.debug.value(x[k])
                    if k < self.N:
                        U_opt[:, k] = opti.debug.value(u[k])
                
                cost = opti.debug.value(opti.f)
                return X_opt, U_opt, cost, False, None
            except:
                return None, None, None, False, None

def main():
    """
    Main function to run trajectory optimization example
    """
    # Create optimizer
    optimizer = CarTrajectoryOptimizer()
    
    print("Car Trajectory Optimization with Spline Reference Tracking")
    print("=" * 60)
    print(f"Time horizon: {optimizer.T:.1f} s")
    print(f"Control intervals: {optimizer.N}")
    print()
    
    # ========== SPLINE REFERENCE TRACKING EXAMPLE ==========
    print("SPLINE REFERENCE TRACKING EXAMPLE")
    print("=" * 60)
    
    # Define waypoints for spline reference trajectory
    waypoints = [
        (0, 0),      # Start
        (5, 3),      # Waypoint 1
        (10, 8),     # Waypoint 2
        (15, 10),    # Waypoint 3
        (20, 12),    # Waypoint 4
        (25, 8),     # Waypoint 5
        # (30, 5),     # Waypoint 6
        # (35, 2),     # Waypoint 7
        # (40, 0),     # End
    ]
    
    desired_speed = 5.0  # m/s
    
    print(f"Waypoints: {waypoints}")
    print(f"Desired speed: {desired_speed:.1f} m/s")
    print()
    
    # Generate spline reference trajectory
    print("Generating spline reference trajectory...")
    reference = optimizer.generate_spline_reference(waypoints, desired_speed)
    print(f"Reference trajectory generated with {len(reference['x'])} points")
    print(f"Total path length: {reference['total_length']:.2f} m")
    print()
    
    # Define start state (close to but not exactly on the reference)
    start_state = np.array([0.5, -0.5, np.pi/6, 1.0])  # slightly off the reference
    
    # Define obstacles that the vehicle must avoid while tracking the reference
    obstacles = [
        (2, 8, 6),   # obstacle near waypoint 7
        (12, 8, 2.0),   # obstacle near waypoint 2
        (22, 10, 1.5),  # obstacle near waypoint 4
    ]
    
    print(f"Start state: x={start_state[0]:.1f}, y={start_state[1]:.1f}, θ={start_state[2]*180/np.pi:.1f}°, v={start_state[3]:.1f} m/s")
    print(f"Number of obstacles: {len(obstacles)}")
    print()
    
    # Solve reference tracking optimization
    print("REFERENCE TRACKING OPTIMIZATION")
    print("-" * 40)
    
    result = optimizer.solve_reference_tracking(start_state, obstacles)
    X_opt, U_opt, cost, success = result[:4]
    
    if X_opt is not None:
        if success:
            print(f"Reference tracking optimization completed successfully!")
        else:
            print(f"Reference tracking optimization did not converge, showing debug solution...")
        
        print(f"Final cost: {cost:.4f}")
        print(f"Final position: x={X_opt[0, -1]:.2f}, y={X_opt[1, -1]:.2f}")
        print(f"Final heading: {X_opt[2, -1]*180/np.pi:.1f}°")
        print(f"Final velocity: {X_opt[3, -1]:.2f} m/s")
        
        # Calculate tracking errors
        pos_errors = np.sqrt((X_opt[0, :] - reference['x'])**2 + (X_opt[1, :] - reference['y'])**2)
        heading_errors = np.abs(X_opt[2, :] - reference['theta']) * 180/np.pi
        heading_errors = np.minimum(heading_errors, 360 - heading_errors)  # Handle wrapping
        velocity_errors = np.abs(X_opt[3, :] - reference['v'])
        
        print(f"Average position error: {np.mean(pos_errors):.3f} m")
        print(f"Max position error: {np.max(pos_errors):.3f} m")
        print(f"Average heading error: {np.mean(heading_errors):.1f}°")
        print(f"Max heading error: {np.max(heading_errors):.1f}°")
        print(f"Average velocity error: {np.mean(velocity_errors):.3f} m/s")
        print(f"Max velocity error: {np.max(velocity_errors):.3f} m/s")
        
        # Plot results
        optimizer.plot_results(X_opt, U_opt, obstacles, success)
        

if __name__ == "__main__":
    main()
