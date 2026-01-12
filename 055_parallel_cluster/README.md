# AWS ParallelCluster with Terraform

Complete Infrastructure as Code solution for AWS ParallelCluster with LAMMPS molecular dynamics simulation.

## Overview

This project provides a fully automated setup for AWS ParallelCluster using Terraform. With a single `terraform apply` command, you can create:
- Complete AWS infrastructure (VPC, subnets, EFS, S3, etc.)
- AWS ParallelCluster with Slurm scheduler
- CPU and GPU compute queues
- Automated software installation
- LAMMPS molecular dynamics framework

## Quick Start (15 minutes)

```bash
# 1. Navigate to terraform directory
cd terraform

# 2. Initialize Terraform
terraform init

# 3. Create everything (infrastructure + cluster)
terraform apply
# Wait ~15 minutes for cluster creation

# 4. Check cluster status
pcluster describe-cluster --cluster-name my-cluster --region ap-northeast-1

# 5. Connect to HeadNode
pcluster ssh --cluster-name my-cluster --region ap-northeast-1

# 6. On HeadNode: Build and test LAMMPS
chmod +x utils/*.sh
./utils/build_lammps_cpu.sh
sbatch utils/run_lammps_mpi_test.sh
```

That's it! Everything is created and ready to use.

## What Gets Created

### Infrastructure Resources

| Resource | Description | Cost |
|----------|-------------|------|
| VPC | 10.0.0.0/16 with DNS hostnames | Free |
| Public Subnet | 10.0.0.0/24 for HeadNode | Free |
| Private Subnet | 10.0.1.0/24 for ComputeNodes | Free |
| Internet Gateway | Public internet access | Free |
| NAT Gateway | Private subnet internet | ~$0.045/hr |
| Route Tables | Public & Private routing | Free |
| S3 VPC Endpoint | Cost-effective S3 access | Free |
| EFS | Shared storage (encrypted) | ~$0.30/GB-month |
| S3 Bucket | Scripts storage | ~$0.025/GB-month |
| SSH Key Pair | ED25519 key (auto-generated) | Free |

### ParallelCluster Resources

| Resource | Description | Cost |
|----------|-------------|------|
| HeadNode | t3.medium (always running) | ~$0.0416/hr |
| CPU Queue | t3.medium x0-2 (auto-scale) | ~$0.0416/hr per node |
| GPU Queue | g4dn.xlarge x0-2 (auto-scale) | ~$0.526/hr per node |

**Total Cost (idle)**: ~$42/month (NAT Gateway + HeadNode only)

### Software Stack

- **OS**: Ubuntu 22.04
- **Scheduler**: Slurm
- **Compiler**: GCC 11.4.0
- **MPI**: OpenMPI 4.1.7
- **Build Tools**: CMake 3.22.1
- **GPU**: CUDA 12.3 (on GPU nodes)
- **Application**: LAMMPS with Kokkos (OpenMP/CUDA)

## Directory Structure

```
055_parallel_cluster/
├── README.md                     # This file
├── config.yaml                   # ParallelCluster config (manual, optional)
├── install_software.sh           # Node initialization script
│
├── terraform/                    # Infrastructure as Code
│   ├── main.tf                   # Provider configuration
│   ├── variables.tf              # Input variables
│   ├── vpc.tf                    # VPC and networking
│   ├── efs.tf                    # EFS file system
│   ├── s3.tf                     # S3 bucket
│   ├── key_pair.tf               # SSH key generation
│   ├── pcluster.tf               # ParallelCluster creation
│   ├── outputs.tf                # Output values
│   ├── config.yaml.tpl           # Config template
│   └── terraform.tfvars.example  # Example variables
│
└── utils/                        # Utility scripts
    ├── build_lammps_cpu.sh       # Build LAMMPS (CPU)
    ├── build_lammps_gpu.sh       # Build LAMMPS (GPU)
    ├── lammps_test_input.lmp    # LAMMPS test input
    ├── run_lammps_mpi_test.sh   # MPI test job
    └── run_lammps_scaling_test.sh # Scaling test
```

## Prerequisites

### Required Software

1. **Terraform** (>= 1.0)
   ```bash
   # macOS
   brew install terraform

   # Linux
   wget https://releases.hashicorp.com/terraform/1.6.0/terraform_1.6.0_linux_amd64.zip
   unzip terraform_1.6.0_linux_amd64.zip
   sudo mv terraform /usr/local/bin/
   ```

2. **AWS CLI** configured with credentials
   ```bash
   aws configure
   # Enter your AWS Access Key ID and Secret Access Key
   ```

3. **AWS ParallelCluster CLI**
   ```bash
   pip3 install aws-parallelcluster
   # Or using uv
   uv pip install aws-parallelcluster
   ```

4. **Verify setup**
   ```bash
   terraform version
   aws sts get-caller-identity
   pcluster version
   ```

## Detailed Setup

### Step 1: Initialize Terraform

```bash
cd terraform
terraform init
```

This downloads required providers: AWS, TLS, Local, Random.

### Step 2: Configure Variables (Optional)

```bash
cp terraform.tfvars.example terraform.tfvars
# Edit terraform.tfvars to customize
```

Default configuration:
- Region: `ap-northeast-1` (Tokyo)
- VPC CIDR: `10.0.0.0/16`
- Cluster name: `my-cluster`

Available options:
```hcl
aws_region = "ap-northeast-1"
cluster_name = "my-hpc-cluster"
```

### Step 3: Plan and Apply

```bash
# Preview changes
terraform plan

# Create everything
terraform apply
# Type 'yes' when prompted

# Expected time:
#   - Infrastructure: 3-5 minutes
#   - ParallelCluster: 10-12 minutes
#   - Total: ~15 minutes
```

### Step 4: Check Cluster Status

Wait for cluster creation to complete (~10-12 minutes):

```bash
# Check cluster status
pcluster describe-cluster \
  --cluster-name my-cluster \
  --region ap-northeast-1

# Or watch status
watch -n 30 "pcluster describe-cluster \
  --cluster-name my-cluster \
  --region ap-northeast-1 \
  --query 'clusterStatus' \
  --output text"
```

Wait until status shows `CREATE_COMPLETE`.

### Step 5: Connect to Cluster

```bash
# Connect to HeadNode via pcluster ssh
pcluster ssh --cluster-name my-cluster --region ap-northeast-1

# This command automatically:
# - Gets the HeadNode IP address
# - Uses the correct SSH key
# - Connects as the ubuntu user
```

Alternatively, you can connect manually if needed:

```bash
# Get HeadNode IP address
pcluster describe-cluster \
  --cluster-name my-cluster \
  --region ap-northeast-1 \
  --query 'headNode.publicIpAddress' \
  --output text

# Connect via SSH manually
ssh -i pcluster-key-ed25519.pem ubuntu@<HeadNode-IP>
```

## LAMMPS Usage

### Build LAMMPS (CPU Version)

On HeadNode:

```bash
# Copy scripts to home directory
cd ~
cp -r /path/to/utils .

# Fix permissions if needed
sudo chown -R ubuntu:ubuntu /shared/lammps

# Build CPU version (~15-20 minutes)
chmod +x utils/build_lammps_cpu.sh
./utils/build_lammps_cpu.sh

# Verify installation
source /shared/lammps/cpu-env.sh
lmp -h
```

### Build LAMMPS (GPU Version)

Submit as Slurm job:

```bash
# Submit GPU build job
sbatch utils/build_lammps_gpu.sh

# Monitor progress
squeue
tail -f /shared/lammps_gpu_build_*.out

# Verify (on GPU node)
source /shared/lammps/gpu-env.sh
lmp -h
```

### Run MPI Parallel Test

```bash
# Create working directory
sudo mkdir -p /shared/lammps_jobs
sudo chown ubuntu:ubuntu /shared/lammps_jobs

# Submit 2-node, 4-task test
sbatch utils/run_lammps_mpi_test.sh

# Monitor
squeue
watch -n 1 sinfo

# View results
cat lammps_test_*.out
```

### Run Scaling Tests

```bash
# Submit scaling test suite
chmod +x utils/run_lammps_scaling_test.sh
./utils/run_lammps_scaling_test.sh

# Tests 4 configurations:
#   - 1 node, 1 task
#   - 1 node, 2 tasks
#   - 2 nodes, 2 tasks
#   - 2 nodes, 4 tasks

# View summary
cat /shared/lammps_scaling_results_*/summary.txt
```

## Cluster Management

### Check Cluster Status

```bash
# Check cluster status
pcluster describe-cluster \
  --cluster-name my-cluster \
  --region ap-northeast-1

# Get specific information
pcluster describe-cluster \
  --cluster-name my-cluster \
  --region ap-northeast-1 \
  --query 'clusterStatus' \
  --output text

# Get HeadNode IP
pcluster describe-cluster \
  --cluster-name my-cluster \
  --region ap-northeast-1 \
  --query 'headNode.publicIpAddress' \
  --output text
```

### View Logs

```bash
# List log streams
pcluster list-cluster-log-streams \
  --cluster-name my-cluster \
  --region ap-northeast-1

# View specific log
pcluster get-cluster-log-events \
  --cluster-name my-cluster \
  --log-stream-name <stream-name> \
  --region ap-northeast-1
```

### Slurm Commands

On HeadNode:

| Command | Description |
|---------|-------------|
| `sinfo` | Show node and partition status |
| `squeue` | Show job queue |
| `sbatch <script>` | Submit batch job |
| `scancel <job_id>` | Cancel job |
| `scontrol show job <job_id>` | Show detailed job info |

## Configuration Options

### Custom Cluster Name

```bash
cat > terraform.tfvars << EOF
cluster_name = "my-hpc-cluster"
EOF

terraform apply
```

### Different Region

```bash
cat > terraform.tfvars << EOF
aws_region = "us-east-1"
availability_zone = "us-east-1a"
EOF

terraform apply
```

### Update Cluster Configuration

1. Edit `terraform/config.yaml.tpl`
2. Run `terraform apply`
3. Cluster will be recreated (EFS data preserved)

## Cleanup

### Delete Everything

```bash
cd terraform
terraform destroy
# Type 'yes' when prompted
```

**What happens:**
1. **ParallelCluster deletion starts** (pcluster delete-cluster)
2. **Wait for deletion** (~5-10 minutes)
   - Terraform automatically waits for cluster deletion to complete
   - You'll see status updates every 10 seconds
3. **Infrastructure deletion** (VPC, EFS, S3, etc.)
4. **Complete cleanup**

**Expected time**: 15-20 minutes total

**Warning**: This deletes EFS data permanently!

### Partial Cleanup

Delete cluster but keep infrastructure:

```bash
pcluster delete-cluster \
  --cluster-name my-cluster \
  --region ap-northeast-1
```

## Troubleshooting

### Cluster Creation Failed

Check logs:
```bash
pcluster list-cluster-log-streams \
  --cluster-name my-cluster \
  --region ap-northeast-1

pcluster get-cluster-log-events \
  --cluster-name my-cluster \
  --log-stream-name <stream> \
  --region ap-northeast-1
```

Common issues:
- VPC DNS hostnames not enabled
- IAM policies missing S3 access
- S3 VPC endpoint not configured
- EFS mount target creation timeout

### SSH Connection Issues

**Recommended: Use pcluster ssh command**

```bash
# This handles everything automatically
pcluster ssh --cluster-name my-cluster --region ap-northeast-1
```

**If using manual SSH connection:**

```bash
# Check cluster status
pcluster describe-cluster --cluster-name my-cluster --region ap-northeast-1

# Verify HeadNode IP
pcluster describe-cluster \
  --cluster-name my-cluster \
  --region ap-northeast-1 \
  --query 'headNode.publicIpAddress' \
  --output text

# Fix key permissions if needed
chmod 400 pcluster-key-ed25519.pem

# Manual SSH
ssh -i pcluster-key-ed25519.pem ubuntu@<HeadNode-IP>
```

### LAMMPS Build Errors

```bash
# Check build logs
cat /shared/lammps_cpu_build.log
cat /shared/lammps_gpu_build.log

# Verify dependencies
gcc --version    # Should be 11.4.0
cmake --version  # Should be 3.22.1
mpirun --version # Should be 4.1.7

# On GPU node
nvcc --version   # Should be 12.3
nvidia-smi       # Check GPU
```

### Permission Denied on /shared

```bash
# Fix ownership
sudo chown -R ubuntu:ubuntu /shared/lammps
sudo chown -R ubuntu:ubuntu /shared/lammps_jobs
```

### Terraform State Issues

```bash
# Refresh state
terraform refresh

# View current state
terraform show

# List resources
terraform state list
```

## Cost Optimization

### Minimize Costs

1. **Use auto-scaling**: Compute nodes scale to 0 when idle
2. **Delete when not in use**: `terraform destroy` removes everything
3. **Choose cheaper regions**: Some regions are less expensive

### Cost Breakdown

**Infrastructure (24/7)**:
- NAT Gateway: $32.40/month (~$0.045/hr)
- Elastic IP (if not attached): $3.60/month

**Compute (when running)**:
- HeadNode (t3.medium): $30/month (~$0.0416/hr, always on)
- CPU node (t3.medium): $0.0416/hr per node (auto-scale)
- GPU node (g4dn.xlarge): $0.526/hr per node (auto-scale)

**Storage**:
- EFS: $0.30/GB-month (only for stored data)
- S3: $0.025/GB-month (negligible for scripts)

**Total (idle with cluster)**: ~$66/month

## Network Architecture

```
┌────────────────────────────────────────────────────────────┐
│ VPC (10.0.0.0/16)                                          │
│                                                            │
│  ┌──────────────────┐  ┌────────────────────┐             │
│  │ Public Subnet    │  │ Private Subnet     │             │
│  │ 10.0.0.0/24      │  │ 10.0.1.0/24        │             │
│  │                  │  │                    │             │
│  │ ┌──────────────┐ │  │ ┌────────────────┐ │             │
│  │ │ HeadNode     │ │  │ │ CPU Queue      │ │             │
│  │ │ t3.medium    │─┼──┼─│ t3.medium x0-2 │ │             │
│  │ │ Ubuntu 22.04 │ │  │ │                │ │             │
│  │ └──────────────┘ │  │ └────────────────┘ │             │
│  │      │ EIP       │  │                    │             │
│  │      │           │  │ ┌────────────────┐ │             │
│  │      │           │  │ │ GPU Queue      │ │             │
│  │      │           │  │ │ g4dn.xlarge    │ │             │
│  │      │           │  │ │ x0-2 (T4 GPU)  │ │             │
│  │      │           │  │ └────────────────┘ │             │
│  │      │           │  │        │           │             │
│  └──────┼───────────┘  └────────┼───────────┘             │
│         │                       │                         │
│         │                  NAT Gateway                    │
│         │                       │                         │
│   Internet Gateway ◄────────────┘                         │
│         │                                                 │
│         │      S3 VPC Endpoint (Gateway)                  │
│         │      EFS (Shared Storage)                       │
└─────────┼─────────────────────────────────────────────────┘
          │
      Internet
```

## Key Features

- ✅ **Infrastructure as Code**: Full Terraform automation
- ✅ **One-Command Setup**: `terraform apply` creates everything
- ✅ **Auto-Scaling**: Compute nodes scale 0-2 based on demand
- ✅ **Dual Queues**: Separate CPU and GPU queues
- ✅ **Ubuntu 22.04**: Modern OS with GCC 11.x
- ✅ **LAMMPS with Kokkos**: Optimized for CPU (OpenMP) and GPU (CUDA)
- ✅ **MPI Support**: Multi-node parallel computing
- ✅ **Shared Storage**: EFS mounted at `/shared`
- ✅ **Cost Optimized**: S3 VPC endpoint, auto-scaling
- ✅ **Reproducible**: Version-controlled infrastructure
- ✅ **Easy Cleanup**: Single `terraform destroy` command

## Advantages Over Manual Setup

### Before (Manual):
```bash
# 20+ AWS CLI commands
aws ec2 create-vpc ...
aws ec2 create-subnet ...
aws ec2 create-internet-gateway ...
# ... many more commands ...
pcluster create-cluster ...
```

### After (Terraform):
```bash
# Single command
terraform apply
```

### Benefits:
1. **Automated**: No manual AWS Console clicking
2. **Reproducible**: Same environment every time
3. **Version Controlled**: Track changes in Git
4. **State Management**: Terraform tracks all resources
5. **Easy Updates**: Modify `.tf` files and apply
6. **Simple Cleanup**: One command to delete everything
7. **Documentation**: Code IS documentation
8. **Team Collaboration**: Share and review changes

## Performance Notes

### CPU Queue (t3.medium)
- 2 vCPUs per instance
- Suitable for testing and small workloads
- LAMMPS with Kokkos OpenMP threading

### GPU Queue (g4dn.xlarge)
- 4 vCPUs, 1x NVIDIA T4 GPU (16GB VRAM)
- CUDA 12.3 with Kokkos CUDA support
- Suitable for GPU-accelerated simulations

### LAMMPS Performance
- **CPU**: Kokkos OpenMP multi-threading
- **GPU**: Kokkos CUDA acceleration
- **MPI**: Multi-node distributed parallelism
- **Hybrid**: MPI + OpenMP/CUDA per node

## Additional Resources

- [Terraform AWS Provider](https://registry.terraform.io/providers/hashicorp/aws/latest/docs)
- [AWS ParallelCluster Documentation](https://docs.aws.amazon.com/parallelcluster/)
- [LAMMPS Documentation](https://docs.lammps.org/)
- [Kokkos Documentation](https://kokkos.github.io/kokkos-core-wiki/)
- [Slurm Documentation](https://slurm.schedmd.com/documentation.html)

## Support and Contribution

For issues and improvements:
1. Check troubleshooting section above
2. Review Terraform state and logs
3. Check AWS ParallelCluster documentation
4. Open an issue with detailed error messages

## License

This project is provided as-is for educational and research purposes.
