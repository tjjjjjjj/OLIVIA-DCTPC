-- Code to set up the 4shooter database for environmental logging.
-- Right now, this is a test database.
-- James Battat
-- February 4, 2010
-- 

USE DMTPC_TEST;

CREATE TABLE IF NOT EXISTS pressure
( 
  	value_cdg FLOAT NOT NULL, 
	rms_cdg FLOAT DEFAULT 0, 
  	value_bpg FLOAT NOT NULL, 
	rms_bpg FLOAT DEFAULT 0, 
  	value_convectron FLOAT NOT NULL, 
	rms_convectron FLOAT DEFAULT 0, 
	setval FLOAT, 
	timestamp TIMESTAMP
);

-- create table wire_hv ( value FLOAT NOT NULL, rms FLOAT DEFAULT 0, setval FLOAT, timestamp TIMESTAMP);
-- create table mesh_hv ( value FLOAT NOT NULL, rms FLOAT DEFAULT 0, setval FLOAT, timestamp TIMESTAMP);
-- create table wire_i ( value FLOAT NOT NULL, rms FLOAT DEFAULT 0, setval FLOAT, timestamp TIMESTAMP);
-- create table mesh_i ( value FLOAT NOT NULL, rms FLOAT DEFAULT 0, setval FLOAT, timestamp TIMESTAMP);
-- create table ccd ( temperature FLOAT NOT NULL, exposure INT NOT NULL, daqtime FLOAT DEFAULT 0, avgpixel FLOAT DEFAULT 0, timestamp TIMESTAMP);
-- create table wire2_hv ( value FLOAT NOT NULL, rms FLOAT DEFAULT 0, setval FLOAT, timestamp TIMESTAMP);
-- create table mesh2_hv ( value FLOAT NOT NULL, rms FLOAT DEFAULT 0, setval FLOAT, timestamp TIMESTAMP);
-- create table wire2_i ( value FLOAT NOT NULL, rms FLOAT DEFAULT 0, setval FLOAT, timestamp TIMESTAMP);
-- create table mesh2_i ( value FLOAT NOT NULL, rms FLOAT DEFAULT 0, setval FLOAT, timestamp TIMESTAMP);
-- create table temp0 ( value FLOAT NOT NULL, rms FLOAT DEFAULT 0, setval FLOAT, timestamp TIMESTAMP);
-- create table temp1 ( value FLOAT NOT NULL, rms FLOAT DEFAULT 0, setval FLOAT, timestamp TIMESTAMP);
-- create table temp2 ( value FLOAT NOT NULL, rms FLOAT DEFAULT 0, setval FLOAT, timestamp TIMESTAMP);
-- alter table ccd add column ccdid INT;
