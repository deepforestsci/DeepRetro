Troubleshooting
===============

.. TODO: Add common issues and their solutions based on user feedback and testing.

This page lists common problems you might encounter and how to resolve them.

Installation Issues
-------------------

- ``ModuleNotFoundError`` for ``aizynthfinder`` (or others):
  - Ensure Conda environment (``dfs_si_challenge``) is activated.
  - Verify ``conda env create -f environment.yml`` completed successfully.
  - Try reinstalling: ``pip install --upgrade --force-reinstall <package_name>``

Configuration / .env File
-------------------------

- Backend ``API_KEY`` not found error:
  - Double-check ``.env`` file is in the project root.
  - Ensure it's named correctly (``.env``).
  - Verify it contains ``API_KEY='yourkey'``.

Frontend Issues
---------------

- Frontend shows errors or can't connect to backend:
  - Ensure backend server (``python src/api.py``) is running.
  - Verify the API key entered in the frontend prompt matches ``.env``.
  - Check browser developer console (F12) for errors.

Model Loading Issues
--------------------

.. TODO: Add potential issues related to model downloading or paths.

Performance Issues
------------------

.. TODO: Add tips for performance issues if applicable.

If your issue is not listed here, please check the GitHub Issues page or consider creating a new one. 