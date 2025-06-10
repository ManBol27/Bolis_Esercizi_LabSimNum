In main1.cpp il codice è parallelizzato ma senza migrazioni (ovvero senza che i core scambino dati tra di loro)

In main2.cpp vengono introdotte anche le migrazioni: presi due core, essi si scambiano l'individuo migliore.