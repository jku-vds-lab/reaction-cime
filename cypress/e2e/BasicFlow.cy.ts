Cypress.on('uncaught:exception', (err, runnable) => {
  return !err.message.includes('scrollLeft');
});

describe('Hydration', () => {
  it('User workflow 1', { scrollBehavior: false }, () => {
    // Upload project
    cy.visit('/');
    cy.createProject('cypress/fixtures/domain_5000_v2.csv');

    // Perform selection on context
    cy.selectSomething();

    // Define groups by clustering
    cy.switchTab('groups');
    cy.byId('define-groups-by-clustering-button').click();
    cy.byId('run-clustering-button').click();

    // Set some encodings
    cy.switchTab('encoding');
    cy.byId('color-encoding-select').click().wait(300).type('measured{downArrow}{enter}');
    cy.byId('size-encoding-select').click();
    cy.chooseNth(3);

    // Open LineUp & Sort after measured_yield
    cy.get('#app').trigger('mousedown', 1000, 1073).wait(300).trigger('mousemove', 1000, 500).wait(300).trigger('mouseup', 1000, 500).wait(1000);
    cy.get('[title="Measured_yield"]').find('[title="Sort"]').filter(':visible').click();

    // Create second view and activate it
    cy.byId('split-view-button').click().wait(1000);
    cy.get('#app').click(1700, 350);

    // Select another encoding for the second view
    cy.switchTab('encoding');
    cy.byId('color-encoding-select').click().wait(300).type('measured{downArrow}{enter}');

    // Set some aggregation for the second view
    cy.byId('custom-tab-1').click();
    cy.byId('aggregate-select').click().wait(300).type('predict{downArrow}{enter}');

    // Activate the second cycle
    cy.get('.MuiSlider-root').filter(':visible').get('.MuiSlider-markLabel[data-index="2"]').filter(':visible').click();
  });
});
