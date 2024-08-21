Usage Guide
=============================

This guide includes in-depth information on function parameters and code documentation for **InPheRNo-ChIP**. 
The guide is divided into three main sections, corresponding to the three scripts of the **InPheRNo-ChIP** framework:

InPheRNo-ChIP Step 1
---------------------

The first step involves processing RNA-seq and ChIP-seq data. 

**Basic Usage**:

.. code-block:: bash

    python InPheC_Step1.py

**Advanced Usage**:

.. automodule:: InPheC_Step1
   :members:
   :undoc-members:
   :show-inheritance:

InPheRNo-ChIP Step 2
---------------------

The second step involves running the probabilistic graphical model (PGM).

**Basic Usage**:

.. code-block:: bash

    python InPheC_Step2.py

**Advanced Usage**:

.. automodule:: InPheC_Step2
   :members:
   :undoc-members:
   :show-inheritance:

InPheRNo-ChIP Step 3
---------------------

The final step combines and processes the output from the previous steps.

**Basic Usage**:

.. code-block:: bash

    python InPheC_Step3.py

**Advanced Usage**:

.. automodule:: InPheC_Step3
   :members:
   :undoc-members:
   :show-inheritance:


FAQs
----

Include a Frequently Asked Questions section to cover common queries:

1. **Question 1**: *can't wait for the 1st question you have for InPheC!*
   .. code-block:: none
   Answer 1:...

2. **Question 2**:
   .. code-block:: none
   Answer 2:...


Support and Contribution
------------------------

We are committed to providing support for **InPheRNo-ChIP** users and actively encourage contributions to the project.

Getting Support
^^^^^^^^^^^^^^^

If you encounter issues or have questions while using **InPheRNo-ChIP**, there are several ways to get support:

- **Email the Author**: Feel free to reach out to the author of the paper associated with **InPheRNo-ChIP**. The author's email can typically be found in the corresponding research paper or on the project's main webpage.
- **GitHub Issues**: For technical problems or bugs, you can open an issue on the GitHub repository. Visit the `InPheRNo-ChIP Issues page <https://github.com/your-username/InPheRNo-ChIP/issues>`_ to submit an issue.

Contributing to InPheRNo-ChIP
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Your contributions to **InPheRNo-ChIP** are highly valued. There are various ways you can contribute to the project:

- **Reporting Bugs**: If you find a bug, please report it by opening an issue on our GitHub repository. Provide as much detail as possible to help us understand and address the issue.
- **Feature Requests**: Have ideas for new features or improvements? We'd love to hear them! Please file a feature request on our `GitHub Issues page <https://github.com/your-username/InPheRNo-ChIP/issues>`_.
- **Code Contributions**: If you're interested in contributing code, feel free to fork the repository and submit a pull request with your changes. For major changes, please open an issue first to discuss what you would like to change.

We appreciate all forms of feedback and contributions to make **InPheRNo-ChIP** better for everyone!

.. note::

   This guide assumes a basic understanding of python and pymc. If you are new to these, consider reading `PyMC documentation <https://www.pymc.io/projects/docs/en/stable/learn/core_notebooks/pymc_overview.html>`_
