import * as React from 'react';

function MainApp() {
  return <div>hello world</div>;
}

describe('Health check for Cypress component test', () => {
  it('should mount MainApp', () => {
    cy.mount(<MainApp />);
    cy.get('div').should('include.text', 'hello world');
  });
});
