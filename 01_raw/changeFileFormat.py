import scipy.io

read_file_path = '00_SCR Optimization/rawdata/soc-wiki-Vote-youtube1k.mtx'
write_file_path = '00_SCR Optimization/rawdata/youtube1k.mat'
# Load sparse matrix from .mtx file
sparse_matrix = scipy.io.mmread(read_file_path)

# Save the matrix to a .mat file
scipy.io.savemat(write_file_path, {'matrix': sparse_matrix})


