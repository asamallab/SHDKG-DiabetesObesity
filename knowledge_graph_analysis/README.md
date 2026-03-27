## Prerequisites

To run the knowledge graph–based analysis, a Neo4j instance is required. You will need the following connection credentials:

- URI  
- Username  
- Password  

---

## Neo4j Setup Options

### 1. Neo4j Community Edition (Local Installation)
You can download and install Neo4j Community Edition (Debian/Ubuntu) from:  
https://go.neo4j.com/download-thanks.html?edition=community&release=2026.02.3&flavour=deb

### 2. Neo4j Web Edition (Cloud Instance)
Alternatively, you can create a cloud-hosted Neo4j instance by searching for Neo4j in your browser, signing in, and launching a new database instance. This will provide the required URI, username, and password.

---

## Project Folder Structure

This directory contains all scripts and datasets required to reproduce the knowledge graph–based analysis workflow.

- **INPUT/** → Contains triplet files and associated metadata used to construct the knowledge graph.  
- **OUTPUT/** → Stores results generated from downstream analyses.  
