Installation Guide
====================================

This guide provides instructions for setting up **InPheRNo-ChIP** using a Conda environment, specifically tailored for Python 3.10.12 and PyMC 5.9.2.

.. contents:: 
   :local:
   :depth: 1

Prerequisites
-------------

Before installing **InPheRNo-ChIP**, ensure that Conda is installed on your system. Conda simplifies package management and deployment. If Conda is not installed, download it from [Anaconda's official website](https://www.anaconda.com/products/distribution).

Creating a Conda Environment
----------------------------

**InPheRNo-ChIP** requires specific versions of libraries, managed using a Conda environment. The project repository includes an `inpheC_0820.yml` file that specifies these dependencies, including Python 3.10.12 and PyMC 5.9.2.

To set up and activate the Conda environment, follow these steps:

1. **Clone the Repository**:
   
   Clone the **InPheRNo-ChIP** repository and navigate to the project directory:

   .. code-block:: bash

      git clone https://github.com/Emad-COMBINE-lab/InPheRNo-ChIP.git
      cd InPheRNo-ChIP

2. **Create the Conda Environment**:
   
   Use the `inpheC_0820.yml` file, which is pre-configured for Python 3.10.12 and PyMC 5.9.2, to create a new Conda environment:

   .. code-block:: bash

      conda env create -f inpheC_0820.yml

   This command will create a new environment named `inpheC` and install all the required dependencies.

3. **Activate the Conda Environment**:

   Activate the newly created environment:

   .. code-block:: bash

      conda activate inpheC


Updating the Environment
------------------------

Should the `inpheC_0820.yml` file be updated, synchronize your environment with:

.. code-block:: bash

   conda env update -f inpheC_0820.yml --prune

Troubleshooting
---------------

If you encounter installation issues:

- **Update Conda**: Ensure your Conda installation is up to date by running `conda update conda`.
- **Check the Environment File**: Verify that the `inpheC_0820.yml` file you are using is the correct one from the **InPheRNo-ChIP** repository.
- **Consult the Issue Tracker**: For specific errors, review the project's issue tracker to see if others have encountered and resolved similar issues.
- **Recreate the Environment**: Sometimes, issues can be resolved by removing the existing Conda environment and recreating it:
  
  .. code-block:: bash

      conda env remove -n inpheC
      conda env create -f inpheC_0820.yml

- **Open a GitHub Issue**: If your issue persists, please open a ticket on our Github. Include detailed error logs and environment information to help us diagnose and address your issue effectively