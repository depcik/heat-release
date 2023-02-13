# heat-release
Open-Source Energy, Entropy, and Exergy 0-D Heat Release Model for Internal Combustion Engines

If you are interested in using this software program, I have set up a Discord channel for conversations, to keep track of bugs, etc. Feel free to e-mail me at depcik at ku dot edu to receive a link to the Discord channel and start a conversation.

Sample input:
* Input Required for HR2022 - this is an Excel spreadsheet that indicates everything necessary as input
* 18.0NM_DME_ULSD_18ESRRaw - this is the raw in-cylinder pressure data including the combustion event
* YanmarL100vMotoringv2014 - this is the raw in-cylinder motoring data for our engine
* SampleBiodiesel - this is a sample speciation of biodiesel

Folder: CHEMKIN III Data - this contains all the CHEMKIN-style curve-fit data. You need it to run the program.

I would highly suggest running the program as-is to reproduce the data in our, hopefully to be published, paper. 

Notes:
2/13/2023 - The initial version has been uploaded. Currently, it works for gaseous-assisted direct injection combustion. The "massenergy.m" file will need to be modified to handle just standard direct injection combustion. This would be the "else" option starting on line 452. It should be straightforward to get it to work, but we have been concentrating on the gas-assisted combustion efforts here at KU. I would be happy to help modify this file to get it to work, just let me know. 
