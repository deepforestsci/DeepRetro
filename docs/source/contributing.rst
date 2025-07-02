Contributing Guide
=================

How to contribute to RecursiveLLM.

Getting Started
--------------

1. **Fork the repository**
2. **Create a feature branch:**

   .. code-block:: bash

      git checkout -b feature/your-feature-name

3. **Make your changes** (see :doc:`development`)
4. **Run tests:**

   .. code-block:: bash

      python -m pytest tests/

5. **Submit a pull request**

Code Standards
--------------

* Follow PEP 8 style guide
* Use type hints for all functions
* Add docstrings (Google style)
* Keep functions under 50 lines
* Write tests for new features

Commit Messages
--------------

Use conventional commit format:

.. code-block:: bash

   feat: add new LLM model support
   fix: resolve API key validation issue
   docs: update installation instructions
   test: add integration tests for cache
   refactor: simplify molecule validation

Pull Request Process
-------------------

1. **Update documentation** if needed
2. **Add tests** for new functionality
3. **Ensure all tests pass**
4. **Update changelog** if applicable
5. **Request review** from maintainers

Issue Reporting
--------------

When reporting issues:

* Use the issue template
* Include error messages and stack traces
* Provide steps to reproduce
* Specify your environment (OS, Python version, etc.)

Testing
-------

Run tests before submitting:

.. code-block:: bash

   # Run all tests
   python -m pytest tests/

   # Run specific test file
   python -m pytest tests/test_api.py

   # Run with coverage
   python -m pytest tests/ --cov=src

Code Review
----------

All PRs require:

* Passing tests
* Code review approval
* Documentation updates
* No merge conflicts

Release Process
--------------

1. **Update version** in `pyproject.toml`
2. **Update changelog**
3. **Create release** on GitHub
4. **Tag release** with version number

Contact
-------

* GitHub Issues: For bugs and feature requests
* Discussions: For questions and ideas
* Email: For security issues 