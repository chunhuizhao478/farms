import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import scipy.special as sp_gamma
from scipy.spatial.distance import cdist
from scipy.interpolate import griddata
import netCDF4
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve
import psutil
import gc

class MemoryEfficientWeibull3D:
    def __init__(self, nx=80, ny=20, nz=10, 
                 x_range=(-4000, 4000), y_range=(-2000, 2000), z_range=(0, 1000),
                 correlation_length=200, chunk_size=1000):
        """
        Initialize 3D Weibull field generator with memory optimization
        
        Parameters:
        -----------
        nx, ny, nz : int
            Grid dimensions
        x_range, y_range, z_range : tuple
            Coordinate ranges
        correlation_length : float
            Spatial correlation length
        chunk_size : int
            Size of chunks for memory-efficient processing
        """
        self.nx, self.ny, self.nz = nx, ny, nz
        self.x_range, self.y_range, self.z_range = x_range, y_range, z_range
        self.correlation_length = correlation_length
        self.chunk_size = chunk_size
        
        # Create coordinate arrays
        self.x = np.linspace(x_range[0], x_range[1], nx)
        self.y = np.linspace(y_range[0], y_range[1], ny)
        self.z = np.linspace(z_range[0], z_range[1], nz)
        
        # Create 3D grid
        self.x_grid, self.y_grid, self.z_grid = np.meshgrid(self.x, self.y, self.z, indexing='ij')
        self.grid_points = np.column_stack((
            self.x_grid.ravel(), 
            self.y_grid.ravel(), 
            self.z_grid.ravel()
        ))
        
        self.num_points = self.grid_points.shape[0]
        print(f"Total grid points: {self.num_points}")
        self._check_memory_requirements()
    
    def _check_memory_requirements(self):
        """Check if we have enough memory for the covariance matrix"""
        matrix_size_gb = (self.num_points ** 2 * 8) / (1024**3)  # 8 bytes per float64
        available_memory_gb = psutil.virtual_memory().available / (1024**3)
        
        print(f"Covariance matrix size: {matrix_size_gb:.2f} GB")
        print(f"Available memory: {available_memory_gb:.2f} GB")
        
        if matrix_size_gb > available_memory_gb * 0.5:
            print("WARNING: Matrix may be too large for available memory!")
            print("Consider reducing grid size or using chunked processing")
    
    def compute_covariance_chunked(self, use_sparse=True, sparsity_threshold=0.01):
        """
        Compute covariance matrix in chunks to save memory
        
        Parameters:
        -----------
        use_sparse : bool
            Whether to use sparse matrix representation
        sparsity_threshold : float
            Threshold below which values are set to zero for sparse matrix
        """
        print("Computing covariance matrix in chunks...")
        
        if use_sparse:
            # For sparse matrix, we'll store only significant correlations
            row_indices = []
            col_indices = []
            data = []
            
            for i in range(0, self.num_points, self.chunk_size):
                end_i = min(i + self.chunk_size, self.num_points)
                chunk_points_i = self.grid_points[i:end_i]
                
                for j in range(0, self.num_points, self.chunk_size):
                    end_j = min(j + self.chunk_size, self.num_points)
                    chunk_points_j = self.grid_points[j:end_j]
                    
                    # Compute distances for this chunk
                    distances = cdist(chunk_points_i, chunk_points_j, metric='euclidean')
                    covariances = np.exp(-distances / self.correlation_length)
                    
                    # Keep only significant correlations
                    mask = covariances > sparsity_threshold
                    if np.any(mask):
                        rows, cols = np.where(mask)
                        row_indices.extend(rows + i)
                        col_indices.extend(cols + j)
                        data.extend(covariances[mask])
                    
                    # Force garbage collection
                    del distances, covariances
                    gc.collect()
            
            # Create sparse matrix
            self.covariance_matrix = csc_matrix(
                (data, (row_indices, col_indices)), 
                shape=(self.num_points, self.num_points)
            )
            print(f"Sparse matrix density: {self.covariance_matrix.nnz / self.num_points**2:.4f}")
            
        else:
            # Dense matrix approach (memory intensive)
            self.covariance_matrix = np.zeros((self.num_points, self.num_points))
            
            for i in range(0, self.num_points, self.chunk_size):
                end_i = min(i + self.chunk_size, self.num_points)
                chunk_points = self.grid_points[i:end_i]
                
                # Compute distances for entire chunk at once
                distances = cdist(chunk_points, self.grid_points, metric='euclidean')
                self.covariance_matrix[i:end_i, :] = np.exp(-distances / self.correlation_length)
                
                del distances
                gc.collect()
    
    def generate_gaussian_field_efficient(self, method='cholesky_chunked'):
        """
        Generate Gaussian field using memory-efficient methods
        
        Parameters:
        -----------
        method : str
            'cholesky_chunked' : Chunked Cholesky decomposition
            'iterative' : Iterative solver (for very large problems)
            'svd_truncated' : Truncated SVD approximation
        """
        print(f"Generating Gaussian field using {method} method...")
        
        if method == 'cholesky_chunked':
            return self._cholesky_chunked()
        elif method == 'iterative':
            return self._iterative_solver()
        elif method == 'svd_truncated':
            return self._svd_truncated()
        else:
            raise ValueError(f"Unknown method: {method}")
    
    def _cholesky_chunked(self):
        """Chunked Cholesky decomposition for memory efficiency"""
        # Add regularization for numerical stability
        if hasattr(self.covariance_matrix, 'toarray'):
            # Sparse matrix
            regularization = 1e-6 * csc_matrix(np.eye(self.num_points))
            regularized_cov = self.covariance_matrix + regularization
        else:
            # Dense matrix
            regularized_cov = self.covariance_matrix + 1e-6 * np.eye(self.num_points)
        
        try:
            # Attempt Cholesky decomposition
            if hasattr(regularized_cov, 'toarray'):
                # Convert to dense for Cholesky (might be memory intensive)
                print("Converting sparse matrix to dense for Cholesky...")
                dense_cov = regularized_cov.toarray()
                L = np.linalg.cholesky(dense_cov)
                del dense_cov
            else:
                L = np.linalg.cholesky(regularized_cov)
            
            # Generate uncorrelated random variables
            z = np.random.randn(self.num_points)
            
            # Generate correlated field
            gaussian_field = L @ z
            
            del L
            gc.collect()
            
            return gaussian_field
            
        except np.linalg.LinAlgError:
            print("Cholesky failed, falling back to SVD method")
            return self._svd_truncated()
    
    def _svd_truncated(self, n_components=None):
        """Truncated SVD approximation for very large matrices"""
        if n_components is None:
            n_components = min(1000, self.num_points // 4)
        
        print(f"Using truncated SVD with {n_components} components...")
        
        if hasattr(self.covariance_matrix, 'toarray'):
            cov_dense = self.covariance_matrix.toarray()
        else:
            cov_dense = self.covariance_matrix
        
        # Compute truncated SVD
        U, s, Vt = np.linalg.svd(cov_dense, full_matrices=False)
        
        # Keep only top components
        U_trunc = U[:, :n_components]
        s_trunc = s[:n_components]
        
        # Generate correlated field using truncated decomposition
        z = np.random.randn(n_components)
        gaussian_field = U_trunc @ (np.sqrt(s_trunc) * z)
        
        del U, s, Vt, U_trunc, s_trunc
        gc.collect()
        
        return gaussian_field
    
    def _iterative_solver(self):
        """Iterative solver for very sparse matrices"""
        print("Using iterative solver...")
        
        # Generate target vector
        target = np.random.randn(self.num_points)
        
        # Solve Cx = target where C is covariance matrix
        # This gives us x such that when we compute Lx (where LL^T = C), we get a correlated field
        if hasattr(self.covariance_matrix, 'toarray'):
            # Sparse solver
            from scipy.sparse.linalg import cg
            gaussian_field, info = cg(self.covariance_matrix, target, maxiter=1000)
            if info != 0:
                print(f"Warning: Iterative solver did not converge (info={info})")
        else:
            # Dense solver
            gaussian_field = np.linalg.solve(self.covariance_matrix, target)
        
        return gaussian_field
    
    def transform_to_weibull(self, gaussian_field, mu=46.8, shape_param=12):
        """Transform Gaussian field to Weibull distribution"""
        print("Transforming to Weibull distribution...")
        
        # Compute scale parameter
        scale_param = mu / sp_gamma.gamma(1 + 1/shape_param)
        
        # Convert to uniform using CDF
        uniform_values = stats.norm.cdf(gaussian_field)
        
        # Apply inverse Weibull CDF
        weibull_field = scale_param * (-np.log(1 - uniform_values))**(1 / shape_param)
        
        return weibull_field, scale_param
    
    # def compute_xi_o(self, weibull_field):
    #     """Compute xi_o transformation"""
    #     print("Computing xi_o transformation...")
        
    #     xi_o_field = -np.sqrt(2) / np.sqrt(
    #         1 + (1 + 1) * (1 + 1) * 
    #         np.sin(weibull_field * np.pi / 180) * 
    #         np.sin(weibull_field * np.pi / 180)
    #     )
        
    #     return xi_o_field
    
    def visualize_3d_slices(self, field, title="3D Field", save_plots=False):
        """Visualize 3D field as 2D slices"""
        field_3d = field.reshape(self.nx, self.ny, self.nz)
        
        # Plot slices at different z levels
        n_slices = min(4, self.nz)
        z_indices = np.linspace(0, self.nz-1, n_slices, dtype=int)
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.ravel()
        
        for i, z_idx in enumerate(z_indices):
            ax = axes[i]
            im = ax.contourf(self.x_grid[:, :, z_idx], self.y_grid[:, :, z_idx], 
                           field_3d[:, :, z_idx], cmap='viridis')
            ax.set_title(f'{title} - Z = {self.z[z_idx]:.1f}')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            plt.colorbar(im, ax=ax)
        
        plt.tight_layout()
        if save_plots:
            plt.savefig(f'{title.lower().replace(" ", "_")}_slices.png', dpi=150, bbox_inches='tight')
        plt.show()
    
    def interpolate_to_netcdf_points(self, field, file_path, output_csv="mapped_weibull_field_3d.csv"):
        """Interpolate field values to NetCDF coordinate points"""
        print("Interpolating to NetCDF points...")
        
        try:
            nc = netCDF4.Dataset(file_path)
            x_coord = nc.variables['coordx'][:]
            y_coord = nc.variables['coordy'][:]
            
            # Check if z coordinates exist
            if 'coordz' in nc.variables:
                z_coord = nc.variables['coordz'][:]
                target_points = np.column_stack((x_coord, y_coord, z_coord))
                
                # Check bounds for 3D
                in_bounds = ((x_coord >= self.x_range[0]) & (x_coord <= self.x_range[1]) & 
                           (y_coord >= self.y_range[0]) & (y_coord <= self.y_range[1]) &
                           (z_coord >= self.z_range[0]) & (z_coord <= self.z_range[1]))
                
                # Interpolate
                mapped_values = np.zeros_like(x_coord)
                if np.any(in_bounds):
                    mapped_values[in_bounds] = griddata(
                        self.grid_points, field, target_points[in_bounds], 
                        method='linear', fill_value=0
                    )
                
                output_data = np.column_stack((x_coord, y_coord, z_coord, mapped_values))
                
            else:
                # Fall back to 2D if no z coordinates
                print("No z coordinates found, using 2D interpolation...")
                target_points = np.column_stack((x_coord, y_coord))
                
                # Use middle z-slice for 2D interpolation
                mid_z_idx = self.nz // 2
                field_2d = field.reshape(self.nx, self.ny, self.nz)[:, :, mid_z_idx].ravel()
                points_2d = self.grid_points[::self.nz][:len(field_2d)]  # Sample every nz points
                
                in_bounds = ((x_coord >= self.x_range[0]) & (x_coord <= self.x_range[1]) & 
                           (y_coord >= self.y_range[0]) & (y_coord <= self.y_range[1]))
                
                mapped_values = np.zeros_like(x_coord)
                if np.any(in_bounds):
                    mapped_values[in_bounds] = griddata(
                        points_2d, field_2d, target_points[in_bounds], 
                        method='linear', fill_value=0
                    )
                
                output_data = np.column_stack((x_coord, y_coord, mapped_values))
            
            # Save results
            np.savetxt(output_csv, output_data, delimiter=",", comments="")
            print(f"Mapped data saved to {output_csv}")
            
            nc.close()
            
        except Exception as e:
            print(f"Error reading NetCDF file: {e}")
            print("Skipping NetCDF interpolation...")
    
    def generate_field(self, file_path='./static_solve_out.e', method='svd_truncated', 
                      mu=46.8, shape_param=12, visualize=True, save_csv=True):
        """
        Complete workflow to generate 3D Weibull field
        
        Parameters:
        -----------
        file_path : str
            Path to NetCDF file
        method : str
            Method for Gaussian field generation
        mu : float
            Mean friction angle
        shape_param : float
            Weibull shape parameter
        visualize : bool
            Whether to create plots
        save_csv : bool
            Whether to save interpolated results
        """
        print("Starting 3D Weibull field generation...")
        
        # Step 1: Compute covariance matrix
        self.compute_covariance_chunked(use_sparse=True)
        
        # Step 2: Generate Gaussian field
        gaussian_field = self.generate_gaussian_field_efficient(method=method)
        
        # Step 3: Transform to Weibull
        weibull_field, scale_param = self.transform_to_weibull(gaussian_field, mu, shape_param)
        
        # Step 4: Compute xi_o
        # xi_o_field = self.compute_xi_o(weibull_field)
        
        print(f'Shape parameter k: {shape_param}')
        print(f'Scale parameter λ: {scale_param:.4f}')
        
        # Step 5: Visualize
        if visualize:
            self.visualize_3d_slices(gaussian_field, "Gaussian Random Field")
            self.visualize_3d_slices(weibull_field, "Weibull Random Field (Friction Angle)")
            self.visualize_3d_slices(xi_o_field, "Xi_o Field")
        
        # Step 6: Interpolate and save
        if save_csv:
            self.interpolate_to_netcdf_points(xi_o_field, file_path)
        
        return gaussian_field, weibull_field, xi_o_field


# Example usage
if __name__ == "__main__":
    # Initialize with reasonable 3D grid size
    generator = MemoryEfficientWeibull3D(
        nx=100, ny=10, nz=10,  # Reduced size for memory efficiency
        correlation_length=200,
        chunk_size=500
    )
    
    # Generate the field
    gaussian_field, weibull_field, xi_o_field = generator.generate_field(
        method='svd_truncated',  # Most memory efficient for large problems
        visualize=True,
        save_csv=True,
        mu=0.1, 
        shape_param=12
    )
    
    print("3D Weibull field generation completed!")