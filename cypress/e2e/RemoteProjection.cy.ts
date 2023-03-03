Cypress.on('uncaught:exception', (err, runnable) => {
  return !err.message.includes('scrollLeft');
});

describe('Remote Projection', () => {
  /**
   * This tests all remote projections.
   */
  it('All Remote Projections', { scrollBehavior: false }, () => {
    cy.visit('/');
    /* TODO: These are disabled as they run way too long in the CI
    cy.createProject('cypress/fixtures/domain_5000_v2.csv');

    cy.switchTab('projection');

    for (const type of ['umapRemote', 'tsneRemote', 'pcaRemote', 'rmOverlap']) {
      cy.byId(`embedding-${type}`).click();
      if (type !== 'rmOverlap') {
        cy.byId('projection-iterations-number-input').click().type('{selectall}').type('5');
      }
      cy.contains('button', 'Start').click();
      cy.contains('div', /5\/5/, { timeout: 60000 }).click();
      cy.byId('projection-control-card-close-button').click();
    }
    */
  });
});
