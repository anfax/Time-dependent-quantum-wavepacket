import os
import paramiko
import getpass

# Define server details
remote_server = '192.168.6.21'
username = 'lih'
password = getpass.getpass('Enter password: ')

# Get the local directory and remote directory from the user
local_directory = os.path.dirname(os.path.abspath(__file__))
remote_base_dir = '/scratch/lih/'
remote_new_dir_name = input('Enter the new remote directory name: ')
remote_directory = os.path.join(remote_base_dir, remote_new_dir_name)

# Create SSH client
ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

try:
    # Attempt to connect to the remote server
    ssh.connect(remote_server, username=username, password=password)
except paramiko.AuthenticationException:
    print("Authentication failed, please verify your credentials")
    exit(1)
except paramiko.SSHException as sshException:
    print(f"Unable to establish SSH connection: {sshException}")
    exit(1)
except Exception as e:
    print(f"Operation error: {e}")
    exit(1)

# Create SFTP client
sftp = ssh.open_sftp()

def upload_directory(local_dir, remote_dir):
    if not os.path.isdir(local_dir):
        raise ValueError(f"{local_dir} is not a directory")
    
    for root, dirs, files in os.walk(local_dir):
        for dirname in dirs:
            local_path = os.path.join(root, dirname)
            relative_path = os.path.relpath(local_path, local_dir)
            remote_path = os.path.join(remote_dir, relative_path)
            try:
                sftp.mkdir(remote_path)
            except IOError:
                pass  # Directory already exists
        
        for filename in files:
            local_path = os.path.join(root, filename)
            relative_path = os.path.relpath(local_path, local_dir)
            remote_path = os.path.join(remote_dir, relative_path)
            sftp.put(local_path, remote_path)

# Ensure the remote base directory exists
try:
    sftp.stat(remote_base_dir)
except IOError:
    print(f"Remote base directory {remote_base_dir} does not exist.")
    exit(1)

# Ensure the remote directory exists
try:
    sftp.mkdir(remote_directory)
except IOError:
    print(f"Remote directory {remote_directory} already exists or could not be created.")
    exit(1)

# Upload the directory
print(f'Send {local_directory}/* -> {remote_directory}')
upload_directory(local_directory, remote_directory)

# Close the SFTP client
sftp.close()

# Submit the job using SLURM
slurm_script = os.path.join(remote_directory, 'sba.sh')
exe = os.path.join(remote_directory, 'a.out')

# Execute commands on the remote server
commands = [
    f'chmod +x {slurm_script}',
    f'chmod +x {exe}',
    f'cd {remote_directory} && sbatch {slurm_script}'
]

for command in commands:
    stdin, stdout, stderr = ssh.exec_command(command)
    print(stdout.read().decode())
    print(stderr.read().decode())

# Close the SSH connection
ssh.close()