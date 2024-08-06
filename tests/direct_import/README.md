# DirectImport for COSMO-SAC

## Introduction
The `DirectImport` feature is a custom enhancement to the [COSMO-SAC package](https://github.com/usnistgov/COSMOSAC), designed to provide more flexibility and organization in handling sigma profiles for different components. Unlike the standard import methods (`VirginiaTechProfileDatabase`, `DelawareProfileDatabase`) in COSMO-SAC, which require all sigma profiles to be stored in one folder and listed in a TXT file, `DirectImport` allows users to store sigma profiles in separate directories. This feature enables specifying the path and name of the sigma file for each component individually, offering a significant improvement in directory structure and accessibility.

## Features
- **Enhanced Organization:** With `DirectImport`, users can separate sigma profiles into different folders, such as one for polymers and another for APIs (Active Pharmaceutical Ingredients), and provide their paths separately. This organization is particularly beneficial for maintaining a clear and structured directory.
- **Flexible Testing:** `DirectImport` facilitates the testing of different sigma profiles for the same component by allowing users to select which sigma profile to use explicitly. This flexibility is invaluable for research and development purposes, where multiple iterations of sigma profiles may need to be evaluated.

