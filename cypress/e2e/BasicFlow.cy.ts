Cypress.on('uncaught:exception', (err, runnable) => {
  return !err.message.includes('scrollLeft');
});

describe('Hydration', () => {
  it('User workflow 1', { scrollBehavior: false }, () => {
    // Upload project
    cy.visit('/');
    cy.createProject('Hydration', 'cypress/fixtures/domain_5000_v2.csv');

    // Test 2 - project UMAP
    /**cy.switchTab('projection');
    cy.byId('embedding-umapRemote').click();
    cy.byId('projection-iterations-number-input').click().type('{selectall}').type('5');
    cy.contains('button', 'Start').click();
    cy.contains('div', /5\/5/, { timeout: 20000 }).click();
**/
    cy.switchTab('encoding');

    // cy.byId('color-encoding-select').click().wait(300).type('measured{downArrow}{enter}');
    // cy.byId('size-encoding-select').click();
    // cy.chooseNth(3);

    cy.get('#app').trigger('mousedown', 1000, 1073).wait(300).trigger('mousemove', 1000, 500).wait(300).trigger('mouseup', 1000, 500).wait(1000);
    cy.get('[title="Measured_yield"]').find('[title="Sort"]').filter(':visible').click();

    // Create second view
    cy.byId('split-view-button').click().wait(1000);

    // Activate it
    cy.get('#app').click(1700, 350);
    cy.switchTab('encoding');
    cy.byId('color-encoding-select').click().wait(300).type('measured{downArrow}{enter}');

    cy.byId('custom-tab-1').click();
    cy.byId('aggregate-select').click().wait(300).type('predict{downArrow}{enter}');

    cy.get('.MuiSlider-root').filter(':visible').get('.MuiSlider-markLabel[data-index="2"]').filter(':visible').click();
  });
});
